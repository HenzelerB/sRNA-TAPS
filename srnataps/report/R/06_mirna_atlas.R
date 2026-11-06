#!/usr/bin/env Rscript
# =============================================================================
# 06_mirna_atlas.R — Figure 7: complete miRNA modification-site atlas
# =============================================================================
suppressPackageStartupMessages(library(optparse))

.srnataps_r_dir <- tryCatch({
  .a <- commandArgs(FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  if (length(.f) > 0) dirname(normalizePath(.f[1])) else getwd()
}, error = function(e) getwd())
source(file.path(Sys.getenv("SRNATAPS_R_DIR", .srnataps_r_dir), "00_setup.R"))

option_list <- list(
  make_option("--outdir", type = "character"),
  make_option("--figdir", type = "character", default = NULL),
  make_option("--gtf", type = "character", default = NULL,
              help = "Ensembl GTF supplied by the workflow")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$outdir)) stop("--outdir is required")
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

replicate_dir <- file.path(opt$outdir, "07d.replicate_calls", "miRNA")
if (!dir.exists(replicate_dir)) {
  legacy_replicate_dir <- file.path(opt$outdir, "07d.replicate_calls_bayesian", "miRNA")
  if (dir.exists(legacy_replicate_dir)) replicate_dir <- legacy_replicate_dir
}
call_files <- list.files(
  replicate_dir, pattern = "^treat_.*_replicate_miRNA_taps\\.tsv$",
  full.names = TRUE
)
gtf_files <- if (!is.null(opt$gtf) && file.exists(opt$gtf)) {
  opt$gtf
} else {
  list.files(file.path(opt$outdir, "04a.genome"),
             pattern = "\\.gtf$", full.names = TRUE)
}

if (length(call_files) == 0) stop("No replicate miRNA call files found")
if (length(gtf_files) == 0) stop("No Ensembl GTF found in 04a.genome")

# Use the system's streaming grep to avoid loading the 1.5 GB GTF into R.
read_mirna_genes <- function(gtf_file) {
  tmp <- tempfile(fileext = ".gtf")
  on.exit(unlink(tmp))
  status <- system2(
    "grep",
    args = c("-F", shQuote('gene_biotype "miRNA"'), shQuote(gtf_file)),
    stdout = tmp
  )
  if (!status %in% c(0L, 1L)) stop("Failed to filter miRNA records from GTF")
  kept <- readLines(tmp, warn = FALSE)
  kept <- kept[grepl("\\tgene\\t", kept)]
  if (length(kept) == 0) stop("No miRNA gene records found in GTF")

  fields <- strsplit(kept, "\\t", fixed = FALSE)
  get_attr <- function(x, key) {
    vapply(x, function(value) {
      match <- regexpr(paste0(key, ' "[^"]+"'), value)
      if (match[[1]] < 0) return(NA_character_)
      hit <- regmatches(value, match)
      sub('"$', "", sub(paste0('^', key, ' "'), "", hit))
    }, character(1), USE.NAMES = FALSE)
  }
  tibble::tibble(
    chrom = vapply(fields, `[[`, character(1), 1),
    gene_start = as.integer(vapply(fields, `[[`, character(1), 4)),
    gene_end = as.integer(vapply(fields, `[[`, character(1), 5)),
    gene_strand = vapply(fields, `[[`, character(1), 7),
    attributes = vapply(fields, `[[`, character(1), 9)
  ) %>%
    dplyr::mutate(
      gene_id = get_attr(attributes, "gene_id"),
      gene_name = get_attr(attributes, "gene_name"),
      gene_name = dplyr::coalesce(gene_name, gene_id),
      mirna_length = gene_end - gene_start + 1L
    ) %>%
    dplyr::select(-attributes)
}

genes <- read_mirna_genes(gtf_files[[1]])

calls <- lapply(call_files, function(f) {
  cell_line <- sub("^treat_(.*)_replicate_miRNA_taps\\.tsv$", "\\1", basename(f))
  readr::read_tsv(f, show_col_types = FALSE,
                   col_types = readr::cols(start = readr::col_integer(),
                                            end = readr::col_integer(),
                                            padj = readr::col_double())) %>%
    dplyr::mutate(cell_line = cell_line, genomic_position = start + 1L)
}) %>%
  dplyr::bind_rows()

# Assign each genomic cytosine to its Ensembl miRNA gene and calculate the
# nucleotide coordinate from the gene's 5' end (strand aware, 1 based).
mapped <- lapply(seq_len(nrow(calls)), function(i) {
  site <- calls[i, ]
  hit <- genes %>%
    dplyr::filter(chrom == site$chrom,
                  gene_start <= site$genomic_position,
                  gene_end >= site$genomic_position) %>%
    dplyr::slice(1)
  if (nrow(hit) == 0) return(NULL)
  dplyr::bind_cols(site, hit %>% dplyr::select(-chrom)) %>%
    dplyr::mutate(
      mirna_position = dplyr::if_else(
        gene_strand == "+",
        genomic_position - gene_start + 1L,
        gene_end - genomic_position + 1L
      )
    )
}) %>%
  dplyr::bind_rows()

if (nrow(mapped) == 0) stop("No replicate-tested calls overlap annotated miRNA genes")

# Collapse identical sites if more than one cell line is present while retaining
# the strongest supported rate and total evidence for the atlas.
atlas <- mapped %>%
  dplyr::group_by(gene_name, gene_id, chrom, gene_start, gene_end, gene_strand,
                  mirna_length, mirna_position, genomic_position, context) %>%
  dplyr::summarise(
    modification_rate = max(mod_rate, na.rm = TRUE),
    coverage = sum(coverage, na.rm = TRUE),
    padj = min(padj, na.rm = TRUE),
    cell_lines = paste(sort(unique(cell_line)), collapse = ","),
    .groups = "drop"
  ) %>%
  dplyr::arrange(gene_name, mirna_position)

# miRBase supplies mature 5p/3p features that are absent from the Ensembl GTF.
mirbase_search_dirs <- unique(c(
  file.path(opt$outdir, "04a.genome"),
  if (!is.null(opt$gtf) && file.exists(opt$gtf)) dirname(opt$gtf) else character()
))
mirbase_files <- unlist(lapply(mirbase_search_dirs, function(d) {
  if (!dir.exists(d)) return(character())
  list.files(d, pattern = "miRBase.*hsa.*\\.gff3$",
             full.names = TRUE, ignore.case = TRUE)
}), use.names = FALSE)
mature_regions <- tibble::tibble(
  gene_name = character(), mature_name = character(),
  mature_arm = character(), mature_from = integer(), mature_to = integer()
)
if (length(mirbase_files) > 0) {
  gff <- readr::read_tsv(
    mirbase_files[[1]], comment = "#", col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) %>%
    dplyr::transmute(
      chrom = sub("^chr", "", X1), feature = X3,
      mature_start = as.integer(X4), mature_end = as.integer(X5),
      strand = X7, attributes = X9,
      mature_name = sub(".*;Name=([^;]+).*", "\\1", attributes)
    ) %>%
    dplyr::filter(feature == "miRNA")

  gene_keys <- atlas %>%
    dplyr::distinct(gene_name, chrom, gene_start, gene_end, gene_strand,
                    mirna_length)
  mature_regions <- lapply(seq_len(nrow(gene_keys)), function(i) {
    gene <- gene_keys[i, ]
    hit <- gff %>%
      dplyr::filter(
        chrom == gene$chrom, strand == gene$gene_strand,
        mature_start <= gene$gene_end, mature_end >= gene$gene_start
      )
    if (nrow(hit) == 0) return(NULL)
    hit %>%
      dplyr::mutate(
        gene_name = gene$gene_name,
        mature_from = if (gene$gene_strand == "+") {
          mature_start - gene$gene_start + 1L
        } else {
          gene$gene_end - mature_end + 1L
        },
        mature_to = if (gene$gene_strand == "+") {
          mature_end - gene$gene_start + 1L
        } else {
          gene$gene_end - mature_start + 1L
        },
        mature_arm = dplyr::case_when(
          grepl("-5p$", mature_name) ~ "Mature 5p",
          grepl("-3p$", mature_name) ~ "Mature 3p",
          TRUE ~ "Mature miRNA"
        )
      ) %>%
      dplyr::filter(mature_to >= 1, mature_from <= gene$mirna_length) %>%
      dplyr::mutate(
        mature_from = pmax(1L, mature_from),
        mature_to = pmin(gene$mirna_length, mature_to)
      )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(gene_name, mature_name, mature_arm, mature_from, mature_to)
}

# Record whether each reported site falls inside an annotated mature product.
atlas <- lapply(seq_len(nrow(atlas)), function(i) {
  site <- atlas[i, ]
  overlap <- mature_regions %>%
    dplyr::filter(gene_name == site$gene_name,
                  mature_from <= site$mirna_position,
                  mature_to >= site$mirna_position)
  site$mature_miRNA <- if (nrow(overlap) > 0) {
    paste(unique(overlap$mature_name), collapse = ",")
  } else {
    NA_character_
  }
  site
}) %>% dplyr::bind_rows()

# Add mature-arm coordinates for display labels while retaining precursor
# coordinates for the horizontal locus/region layout. Sites outside an
# annotated mature product are intentionally left unlabeled rather than
# displaying misleading precursor positions.
atlas <- lapply(seq_len(nrow(atlas)), function(i) {
  site <- atlas[i, ]
  overlap <- mature_regions %>%
    dplyr::filter(gene_name == site$gene_name,
                  mature_from <= site$mirna_position,
                  mature_to >= site$mirna_position)
  site$mature_position <- if (nrow(overlap) > 0) {
    paste(unique(site$mirna_position - overlap$mature_from + 1L), collapse = ",")
  } else {
    NA_character_
  }
  site
}) %>% dplyr::bind_rows()

readr::write_tsv(atlas, file.path(opt$figdir, "07_mirna_modification_atlas.tsv"))
if (nrow(mature_regions) > 0) {
  readr::write_tsv(
    mature_regions,
    file.path(opt$figdir, "07_mirna_mature_regions_mirbase.tsv")
  )
}

atlas <- atlas %>%
  dplyr::mutate(neg_log10_padj = pmin(-log10(pmax(padj, 1e-300)), 20))

gene_order <- atlas %>%
  dplyr::distinct(gene_name) %>%
  dplyr::arrange(gene_name) %>%
  dplyr::pull(gene_name)
atlas <- atlas %>%
  dplyr::mutate(
    gene_name = factor(gene_name, levels = rev(gene_order))
  )
gene_ranges <- atlas %>%
  dplyr::distinct(gene_name, mirna_length)
mature_regions <- mature_regions %>%
  dplyr::mutate(gene_name = factor(gene_name, levels = rev(gene_order))) %>%
  dplyr::filter(!is.na(gene_name))

p_atlas <- ggplot(atlas, aes(y = gene_name)) +
  geom_segment(
    data = gene_ranges,
    aes(x = 1, xend = mirna_length, yend = gene_name),
    colour = "black", linewidth = 0.5, lineend = "round",
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  geom_segment(
    data = mature_regions,
    aes(x = mature_from, xend = mature_to, y = gene_name, yend = gene_name,
        colour = mature_arm),
    linewidth = 2.2, alpha = 0.5, lineend = "round",
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(x = mirna_position, fill = neg_log10_padj, size = coverage),
    shape = 21, colour = "transparent", alpha = 0.7, stroke = 0
  ) +
  geom_text(
    aes(x = mirna_position, label = mature_position),
    colour = "#1C4062", size = 1.85, vjust = 0.5,
    fontface = "bold",
    check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Mature 5p" = AURORA_MINT,
               "Mature 3p" = AURORA_BLUE,
               "Mature miRNA" = AURORA_GOLD),
    name = "miRBase v22 region"
  ) +
  scale_fill_gradientn(
    colours = c(AURORA_NAVY, AURORA_BLUE, AURORA_GOLD),
    limits = c(0, 20), oob = scales::squish,
    name = expression(-log[10](p[adj])~"(capped at 20)")
  ) +
  scale_size_continuous(
    trans = "sqrt", range = c(2, 6),
    breaks = scales::breaks_pretty(n = 4),
    labels = scales::label_comma(), name = "Coverage"
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    title = "Landscape of miRNA modification sites",
    subtitle = paste0(
      nlevels(atlas$gene_name), " miRNAs | ", nrow(atlas),
      " replicate-tested sites, coloured by significance vs. untreated control"
    ),
    x = "Position within annotated miRNA precursor (nt, 5′ → 3′); labels show mature-arm positions",
    y = NULL,
    caption = paste0(
      "Arrows represent Ensembl GRCh38.112 precursor loci; numbers identify mature-arm positions. ",
      "Sites without a mature-arm annotation are left unlabeled; coordinates are 1-based and strand-normalized."
    )
  ) +
  theme_srnataps(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E6EDF3", linewidth = 0.3),
    axis.text.y = element_text(size = 6.2, face = "italic", colour = "#1C4062"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 2)
  ) +
  guides(
    colour = guide_legend(title.position = "top", title.hjust = 0.5,
                          nrow = 1, order = 1),
    fill = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(38, "mm"), barheight = grid::unit(3, "mm"),
      order = 2
    ),
    size = guide_legend(title.position = "top", title.hjust = 0.5,
                        nrow = 1, order = 3,
                        override.aes = list(fill = AURORA_BLUE,
                                            colour = "transparent", alpha = 0.7))
  )

save_figure(
  p_atlas, file.path(opt$figdir, "07_mirna_modification_atlas.pdf"),
  width = 10, height = max(9, 2.8 + 0.22 * nlevels(atlas$gene_name))
)

message("Figure 7: miRNA modification atlas — ",
        nlevels(atlas$gene_name), " miRNAs, ", nrow(atlas), " sites, ",
        dplyr::n_distinct(mature_regions$gene_name),
        " precursors with miRBase mature annotation — done")
