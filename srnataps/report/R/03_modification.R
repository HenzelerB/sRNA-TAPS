# =============================================================================
# 03_modification.R — Figure 4: detailed TAPS modification landscape
#
# Inputs:
#   - 07.taps_calls/<biotype>/<sample>_<biotype>_taps.tsv
#
# Outputs:
#   - report/figures/04a_modrate_distribution.{pdf,png,svg}
#   - report/figures/04b_top_sites.{pdf,png,svg}
#   - report/figures/04c_condition_comparison.{pdf,png,svg}
#   - report/figures/04d_waterfall.{pdf,png,svg}
#   - report/figures/04e_trinucleotide_context.{pdf,png,svg}
#   - report/figures/04f_metagene_position.{pdf,png,svg}
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
  library(ggridges)
})
# Locate 00_setup.R: env var SRNATAPS_R_DIR wins (set by the pipeline); otherwise
# fall back to this script's own directory so manual runs work on any machine.
.srnataps_r_dir <- tryCatch({
  .a <- commandArgs(FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  if (length(.f) > 0) dirname(normalizePath(.f[1])) else getwd()
}, error = function(e) getwd())
source(file.path(Sys.getenv("SRNATAPS_R_DIR", .srnataps_r_dir), "00_setup.R"))

option_list <- list(
  make_option("--outdir",   type = "character"),
  make_option("--figdir",   type = "character", default = NULL),
  make_option("--min-cov",  type = "integer",   default = 5,
              help = "Minimum coverage to include a site [default: 5]"),
  make_option("--biotypes", type = "character",
              default = "miRNA,tRNA,rRNA,snoRNA",
              help = "Comma-separated biotypes to include in modification plots")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

BIOTYPES_PLOT <- strsplit(opt$biotypes, ",")[[1]]
# The configured pipeline output is 07.taps_calls.  Keep a compatibility
# fallback for older runs that used the former 07.taps_calls_3letter name.
CALLS_DIR     <- file.path(opt$outdir, "07.taps_calls")
if (!dir.exists(CALLS_DIR)) {
  legacy_calls_dir <- file.path(opt$outdir, "07.taps_calls_3letter")
  if (dir.exists(legacy_calls_dir)) CALLS_DIR <- legacy_calls_dir
}

# ── Load all TAPS calls ───────────────────────────────────────────────────────
load_calls <- function(calls_dir, biotypes, min_cov) {
  all_files <- unlist(lapply(biotypes, function(bt) {
    list.files(file.path(calls_dir, bt), pattern = "_taps_annotated\\.tsv$",
               full.names = TRUE)
  }))

  if (length(all_files) == 0) stop("No TAPS call files found in ", calls_dir)
  message("  Loading ", length(all_files), " call files...")

  lapply(all_files, function(f) {
    bt     <- basename(dirname(f))
    sample <- sub(paste0("_", bt, "_taps_annotated\\.tsv$"), "", basename(f))
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE, col_types = cols(chrom = col_character(), start = col_integer(), end = col_integer(), mod_count = col_double(), unmod_count = col_double(), coverage = col_double(), mod_rate = col_double(), pvalue = col_double(), padj = col_double(), snp_flag = col_character(), gene_name = col_character(), gene_id = col_character(), gene_biotype = col_character())) %>%
        dplyr::filter(coverage >= min_cov, snp_flag == "PASS") %>%
        dplyr::mutate(
          biotype   = bt,
          sample    = sample,
          condition = get_condition(sample),
          cell_line = get_cell_line(sample),
          site_label = dplyr::if_else(!is.na(gene_name) & gene_name != ".", gene_name, paste0(chrom, ":", start, "-", end))
        )
      df
    }, error = function(e) NULL)
  }) %>%
    dplyr::bind_rows(Filter(Negate(is.null), .)) %>%
    dplyr::mutate(
      condition = factor(condition, levels = names(CONDITION_COLOURS)),
      biotype   = factor(biotype,   levels = names(BIOTYPE_COLOURS))
    )
}

calls <- load_calls(CALLS_DIR, BIOTYPES_PLOT, opt$`min-cov`)
message("  Total PASS sites loaded: ", nrow(calls))

# ── RNA-relative ("mature") position, from gene span in the GTF ──────────────
# Site labels should read as a position on the transcript (e.g. "pos 34"),
# not a genomic coordinate — chrom:start means nothing to a reader without
# the reference build open. These genes are effectively single-exon at this
# scale (miRNA hairpin / tRNA / snoRNA / small rRNA pseudogene), so gene
# start/end from the GTF stands in for the mature transcript span; position 1
# is the annotated 5' end regardless of genomic strand.
GENOME_DIR <- file.path(dirname(CALLS_DIR), "04a.genome")
GTF_PATH   <- file.path(GENOME_DIR, "Homo_sapiens.GRCh38.112.gtf")

load_gene_spans <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    message("  WARNING: GTF not found at ", gtf_path, " — mature positions unavailable")
    return(dplyr::tibble(gene_id = character(), gene_start = integer(),
                         gene_end = integer(), gene_strand = character()))
  }
  tmp <- tempfile(fileext = ".gtf")
  system2("awk", c("-F", "'\\t'", "'$3==\"gene\"'", shQuote(gtf_path)), stdout = tmp)
  genes <- readr::read_tsv(tmp, col_names = FALSE, show_col_types = FALSE,
                           col_types = readr::cols(.default = "c"))
  unlink(tmp)
  genes %>%
    dplyr::transmute(
      gene_id     = stringr::str_match(X9, 'gene_id "([^"]+)"')[, 2],
      gene_start  = as.integer(X4),
      gene_end    = as.integer(X5),
      gene_strand = X7
    ) %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

gene_spans <- load_gene_spans(GTF_PATH)
message("  Gene spans loaded: ", nrow(gene_spans))

# ── tRNA-specific: Sprinzl/universal tRNA numbering, not linear mature_pos ──
# A simple linear count from a tRNA gene's own annotated start (mature_pos,
# above) is NOT the numbering used in the tRNA modification literature --
# tRNAs share one conserved cloverleaf structure but their loops vary in
# length, so the field aligns every tRNA to a common structural template and
# numbers the SLOTS (Sprinzl numbering), skipping slots a given tRNA lacks
# and adding lettered sub-positions (17A, 27a, ...) where it has an extra
# base. Confirmed 2026-09-07 this matters in practice: this project's own
# linear "mature_pos 52" for MT-TS2 is a real C, but a literature reader
# checking "Sprinzl position 52" finds a totally different base (T/U) at a
# different genomic position (MT:12243, not MT:12258) -- the two numbers
# only coincide by digit, not by nucleotide.
# Reference table built 2026-09-07 from MAMIT-tRNA's own human mt-tRNA
# numbering-conversion tables (Anderson et al. 1981 numbering; mamit-trna.
# u-strasbg.fr, via Wayback Machine archive since the live site is down),
# covering all 22 human mitochondrial tRNA genes.
SPRINZL_PATH <- file.path(dirname(CALLS_DIR), "report", "human_mt_tRNA_sprinzl_numbering.tsv")
load_sprinzl_numbering <- function(path) {
  if (!file.exists(path)) {
    message("  WARNING: Sprinzl numbering table not found at ", path, " — tRNA positions will show linear numbering only")
    return(dplyr::tibble(gene_name = character(), genomic_position_1based = integer(), universal_position = character()))
  }
  readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::mutate(genomic_position_1based = as.integer(genomic_position_1based))
}
sprinzl <- load_sprinzl_numbering(SPRINZL_PATH)
message("  Sprinzl/universal tRNA numbering loaded: ", nrow(sprinzl), " positions across ", dplyr::n_distinct(sprinzl$gene_name), " mt-tRNA genes")

calls <- calls %>%
  dplyr::left_join(gene_spans, by = "gene_id") %>%
  dplyr::left_join(
    dplyr::rename(sprinzl, sprinzl_position = universal_position),
    by = c("gene_name" = "gene_name", "end" = "genomic_position_1based")
  ) %>%
  dplyr::mutate(
    # `start` is 0-based (BED-style, from the calling pipeline); gene_start/
    # gene_end are 1-based (straight from the GTF). `end` (= start + 1 for a
    # single-base feature) IS the true 1-based genomic position, so it -- not
    # `start` -- is what belongs against 1-based annotation coordinates.
    # Verified 2026-09-07: fixed a systematic off-by-one (confirmed by direct
    # FASTA lookup and by disagreeing with the independently-validated MT-TS2
    # universal-position mapping in manuscript_figure_red_thread.md, which
    # places the C48-50 sites at mature positions 33-35, not the 32-34 the
    # unfixed formula gave).
    mature_pos = dplyr::case_when(
      is.na(gene_start)     ~ NA_integer_,
      gene_strand == "+"    ~ end - gene_start + 1L,
      gene_strand == "-"    ~ gene_end - end + 1L,
      TRUE                  ~ NA_integer_
    )
  )

# ── miRNA-specific: true MATURE arm position (5p/3p), not hairpin position ───
# The GTF gene span above is the pri-miRNA hairpin (~60-100nt); the actual
# mature miRNA is a ~20-23nt product cut from one arm of it. miRBase's GFF3
# annotates each mature arm as its own "miRNA" feature (distinct from the
# "miRNA_primary_transcript" hairpin), so a site's mature position + arm
# comes from genomic overlap with that specific short interval, not the
# hairpin span.
MIRBASE_GFF3 <- file.path(GENOME_DIR, "miRBase_v22_hsa_GRCh38.gff3")

load_mirna_mature_arms <- function(gff3_path) {
  if (!file.exists(gff3_path)) {
    message("  WARNING: miRBase GFF3 not found at ", gff3_path, " — mature arm positions unavailable")
    return(dplyr::tibble(chrom = character(), mat_start = integer(), mat_end = integer(),
                         mat_strand = character(), arm = character(), mirna_name = character()))
  }
  lines <- readr::read_lines(gff3_path)
  lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
  gff <- readr::read_tsv(I(paste(lines, collapse = "\n")), col_names = FALSE,
                         show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  gff %>%
    dplyr::filter(X3 == "miRNA") %>%
    dplyr::transmute(
      chrom      = sub("^chr", "", X1),
      mat_start  = as.integer(X4),
      mat_end    = as.integer(X5),
      mat_strand = X7,
      mirna_name = stringr::str_match(X9, "Name=([^;]+)")[, 2],
      arm        = dplyr::case_when(
        grepl("-5p$", mirna_name) ~ "5p",
        grepl("-3p$", mirna_name) ~ "3p",
        TRUE                      ~ NA_character_
      )
    )
}

mirna_arms <- load_mirna_mature_arms(MIRBASE_GFF3)
message("  miRNA mature arms loaded: ", nrow(mirna_arms))

mirna_arm_hits <- calls %>%
  dplyr::filter(biotype == "miRNA") %>%
  dplyr::distinct(chrom, start, end) %>%
  # Same 0-based `start` vs 1-based annotation mismatch as mature_pos above:
  # joining on raw `start` against 1-based mat_start/mat_end wrongly excluded
  # a site sitting exactly on an arm's first base, and wrongly included one
  # extra base past the true 3' end. `end` is the true 1-based position.
  dplyr::inner_join(
    mirna_arms,
    by = dplyr::join_by(chrom == chrom, end >= mat_start, end <= mat_end)
  ) %>%
  dplyr::mutate(
    arm_pos = dplyr::if_else(mat_strand == "+", end - mat_start + 1L, mat_end - end + 1L),
    arm_len = mat_end - mat_start + 1L
  ) %>%
  # A site can in principle sit in >1 overlapping annotated arm; keep one.
  dplyr::distinct(chrom, start, .keep_all = TRUE) %>%
  dplyr::select(chrom, start, arm, arm_pos, arm_len)

calls <- calls %>% dplyr::left_join(mirna_arm_hits, by = c("chrom", "start"))


# ── Figure 3a: Modification rate distribution ────────────────────────────────

calls_dist <- calls %>% dplyr::filter(condition %in% CONDITIONS)

# Histogram — field standard for m5C stoichiometry distributions (binwidth=5%, log10 y)
p_hist <- ggplot(calls_dist,
                 aes(x = mod_rate, fill = condition, colour = condition)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 position = "identity", alpha = 0.6, colour = NA) +
  facet_grid(cell_line ~ biotype) +
  scale_fill_manual(values = CONDITION_COLOURS, labels = CONDITION_LABELS,
                    name = "Condition") +
  scale_colour_manual(values = CONDITION_COLOURS, labels = CONDITION_LABELS,
                      name = "Condition") +
  scale_x_continuous(labels = function(x) x * 100,
                     breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(trans  = "log10",
                     breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000),
                     labels = label_comma(),
                     minor_breaks = NULL) +
  labs(
    title    = "TAPS m5C modification rate distribution",
    subtitle = "Per-site mod_rate at PASS positions (SNP-filtered, BH-corrected)",
    x        = "Modification rate (%)",
    y        = "Number of sites (log10)",
    caption  = paste0("Bin width: 5% | Min coverage: ", opt$`min-cov`,
                      "x | SNP_flag == PASS only")
  ) +
  theme_srnataps() +
  theme(
    text          = element_text(size = 7, family = "Arial", colour = "black"),
    axis.text     = element_text(size = 7, family = "Arial", colour = "black"),
    axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
    plot.title    = element_text(size = 7, family = "Arial", colour = "black"),
    plot.subtitle = element_text(size = 6, family = "Arial", colour = "grey50"),
    legend.text   = element_text(size = 7, family = "Arial", colour = "black"),
    legend.title  = element_text(size = 7, family = "Arial", colour = "black"),
    strip.text    = element_text(size = 7, family = "Arial", colour = "white", face = "bold"),
    plot.caption  = element_text(size = 6, family = "Arial", colour = "grey50")
  )

save_figure(p_hist, file.path(opt$figdir, "04a_modrate_distribution.pdf"),
            width = 6.3, height = 3.6)

message("Figure 4a: Modification rate distribution — done")

# ── Figure 3b: High-confidence site yield per biotype ─────────────────────────
# Not "here are some named genes" (that's Figure 4-6's job for miRNA
# specifically) but "how many confident sites does the assay find per
# biotype" — an aggregate/statistical claim consistent with what 3A/3B/3C/3E
# are each making. Same high-confidence definition used for the sequence-
# context panels (mod_rate >= 50%, padj < 0.01), so the "n = ..." language is
# consistent across the whole figure.
suppressPackageStartupMessages(library(VennDiagram))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_sites <- calls %>%
  dplyr::filter(condition == TREAT_COND, mod_rate >= 0.20, pvalue < 0.05,
                biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA")) %>%
  dplyr::mutate(site_id = paste0(chrom, ":", start)) %>%
  dplyr::distinct(biotype, cell_line, site_id)

BIOTYPES_4 <- c("miRNA", "tRNA", "rRNA", "snoRNA")

if (nrow(venn_sites) > 0 && all(c("Caco2", "HEK") %in% CELL_LINES)) {
  make_venn_panel <- function(bt) {
    caco2 <- venn_sites %>% dplyr::filter(biotype == bt, cell_line == "Caco2") %>% dplyr::pull(site_id)
    hek   <- venn_sites %>% dplyr::filter(biotype == bt, cell_line == "HEK")   %>% dplyr::pull(site_id)
    n_shared <- length(intersect(caco2, hek))

    venn_grob <- grid::grid.grabExpr({
      grid::grid.draw(VennDiagram::draw.pairwise.venn(
        area1 = length(caco2), area2 = length(hek), cross.area = n_shared,
        category = c(paste0("Caco2 (n=", scales::comma(length(caco2)), ")"),
                     paste0("HEK (n=", scales::comma(length(hek)), ")")),
        fill = c(CELL_COLOURS[["Caco2"]], CELL_COLOURS[["HEK"]]),
        alpha = 0.6, lwd = 1, col = "grey40",
        cex = 0.7, cat.cex = 0.55, cat.fontfamily = "Arial", fontfamily = "Arial",
        cat.dist = 0.01, ind = FALSE
        # scaled/euler.d default to TRUE: circle radius is area-proportional
        # to Caco2 vs HEK *within this biotype's own panel* (draw.pairwise.venn
        # always uses a=b, i.e. a true circle, per side — the earlier
        # "stretched ellipse" look was the non-square viewport, now fixed
        # below, not a property of proportional scaling itself). Panels are
        # never compared to each other's absolute scale, only within.
      ))
    })

    gridExtra::arrangeGrob(
      grid::textGrob(bt, gp = grid::gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial")),
      venn_grob,
      ncol = 1, heights = grid::unit(c(1.2, 5), "cm")
    )
  }

  venn_panels <- lapply(BIOTYPES_4, make_venn_panel)

  # Spacer columns on both sides of and between the two panel columns —
  # without them the category labels ("Caco2 (n=...)") of adjacent panels
  # collide, and the leftmost/rightmost labels clip past the canvas edge.
  spacer <- function() grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
  venn_figure <- gridExtra::arrangeGrob(
    grobs = list(
      spacer(), venn_panels[[1]], spacer(), venn_panels[[2]], spacer(),
      spacer(), venn_panels[[3]], spacer(), venn_panels[[4]], spacer()
    ),
    ncol = 5, nrow = 2,
    # Fixed cm widths, not proportional "null" units — the venn cell must be
    # exactly as wide as it is tall (5cm, matching the row height above) for
    # the circles to render as true circles instead of stretched ellipses.
    widths = grid::unit(c(1.8, 5, 1.8, 5, 1.8), "cm"),
    top = grid::textGrob(
      paste0("High-confidence site overlap, Caco2 vs HEK (mod_rate ≥ 20%, p < 0.05, ",
            CONDITION_LABELS[TREAT_COND], ")"),
      gp = grid::gpar(fontsize = 7, fontfamily = "Arial")
    )
  )

  base_out <- file.path(opt$figdir, "04b_top_sites")
  ggplot2::ggsave(paste0(base_out, ".pdf"), venn_figure, width = 6.1, height = 5.2,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(base_out, ".png"), venn_figure, width = 6.1, height = 5.2, dpi = 300)
  ggplot2::ggsave(paste0(base_out, ".svg"), venn_figure, width = 6.1, height = 5.2, device = "svg")
  message("Figure 4b: High-confidence site overlap (Venn, Caco2 vs HEK) — done")

  # ── Figure 3h: Combined Venn — all 4 biotypes pooled, Caco2 vs HEK ─────────
  # Same high-confidence definition and site universe as 3B, but instead of
  # 4 separate per-biotype panels, every miRNA/tRNA/rRNA/snoRNA site is
  # pooled into ONE big Caco2-vs-HEK comparison — the single-number version
  # of "is discovery reproducible across cell lines" rather than the
  # per-biotype breakdown. site_id is biotype-qualified (not just
  # chrom:start) so a coordinate collision across biotypes is never
  # silently merged as "the same site".
  venn_sites_pooled <- venn_sites %>%
    dplyr::mutate(site_id = paste0(biotype, ":", site_id)) %>%
    dplyr::distinct(cell_line, site_id)

  caco2_pooled <- venn_sites_pooled %>% dplyr::filter(cell_line == "Caco2") %>% dplyr::pull(site_id)
  hek_pooled   <- venn_sites_pooled %>% dplyr::filter(cell_line == "HEK")   %>% dplyr::pull(site_id)
  n_shared_pooled <- length(intersect(caco2_pooled, hek_pooled))

  venn_grob_pooled <- grid::grid.grabExpr({
    grid::grid.draw(VennDiagram::draw.pairwise.venn(
      area1 = length(caco2_pooled), area2 = length(hek_pooled), cross.area = n_shared_pooled,
      category = c(paste0("Caco2 (n=", scales::comma(length(caco2_pooled)), ")"),
                   paste0("HEK (n=", scales::comma(length(hek_pooled)), ")")),
      fill = c(CELL_COLOURS[["Caco2"]], CELL_COLOURS[["HEK"]]),
      alpha = 0.6, lwd = 1.2, col = "grey40",
      cex = 1, cat.cex = 0.8, cat.fontfamily = "Arial", fontfamily = "Arial",
      cat.dist = 0.03, ind = FALSE
    ))
  })

  venn_figure_pooled <- gridExtra::arrangeGrob(
    grid::rectGrob(gp = grid::gpar(col = NA, fill = NA)),
    venn_grob_pooled,
    grid::rectGrob(gp = grid::gpar(col = NA, fill = NA)),
    ncol = 3,
    # Same square-viewport trick as 3B (panel must be exactly as wide as
    # tall for draw.pairwise.venn's circles to render as true circles).
    widths = grid::unit(c(1.5, 9, 1.5), "cm"),
    heights = grid::unit(9, "cm"),
    top = grid::textGrob(
      paste0("High-confidence site overlap, Caco2 vs HEK — all biotypes pooled\n(miRNA + tRNA + rRNA + snoRNA, mod_rate ≥ 20%, p < 0.05, ",
            CONDITION_LABELS[TREAT_COND], ")"),
      gp = grid::gpar(fontsize = 9, fontfamily = "Arial")
    )
  )

  base_out_pooled <- file.path(opt$figdir, "04h_top_sites_pooled_venn")
  ggplot2::ggsave(paste0(base_out_pooled, ".pdf"), venn_figure_pooled, width = 5.3, height = 4.6,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(base_out_pooled, ".png"), venn_figure_pooled, width = 5.3, height = 4.6, dpi = 300)
  ggplot2::ggsave(paste0(base_out_pooled, ".svg"), venn_figure_pooled, width = 5.3, height = 4.6, device = "svg")
  message("Figure 4h: High-confidence site overlap (pooled Venn, Caco2 vs HEK, all biotypes) — done")

  # ── Figure 3i: Cell-line concordance — quantitative agreement, not just overlap ─
  # The Venn (3H) shows whether a site clears the HC threshold in BOTH cell
  # lines (binary). This shows whether the two cell lines agree on HOW
  # modified a site is (quantitative) — a site can look "shared" while the
  # two cell lines still disagree substantially on magnitude, and this is the
  # panel that would catch that. No mod_rate floor here (unlike the Venn's
  # >=20% HC threshold): every significant (p<0.05) site present in both
  # cell lines is plotted, so low-rate sites are visible too.
  concordance_sites <- calls %>%
    dplyr::filter(condition == TREAT_COND, pvalue < 0.05,
                  biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA")) %>%
    dplyr::mutate(site_id = paste0(biotype, ":", chrom, ":", start)) %>%
    dplyr::group_by(site_id, biotype, cell_line) %>%
    dplyr::summarise(mod_rate = mean(mod_rate), .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = c(site_id, biotype), names_from = cell_line, values_from = mod_rate) %>%
    dplyr::filter(!is.na(Caco2), !is.na(HEK)) %>%
    dplyr::mutate(biotype = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")))

  if (nrow(concordance_sites) > 0) {
    overall_r    <- stats::cor(concordance_sites$Caco2, concordance_sites$HEK, method = "pearson")
    overall_rmse <- sqrt(mean((concordance_sites$Caco2 - concordance_sites$HEK)^2)) * 100

    p_concordance <- ggplot(concordance_sites, aes(x = Caco2, y = HEK, colour = biotype)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      geom_point(size = 1.6, alpha = 0.55, stroke = 0) +
      annotate("text", x = 0.02, y = 0.98,
              label = paste0("n=", scales::comma(nrow(concordance_sites)),
                            " sites, r=", sprintf("%.2f", overall_r),
                            ", RMSE=", sprintf("%.1f", overall_rmse), "pp"),
              hjust = 0, vjust = 1, size = 1.6, family = "Arial", colour = "grey30") +
      scale_colour_manual(values = BIOTYPE_COLOURS, name = NULL, guide = "none") +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                         breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                         breaks = seq(0, 1, 0.2)) +
      labs(x = "Caco2 mod_rate", y = "HEK mod_rate") +
      theme_srnataps() +
      theme(
        text          = element_text(size = 6, family = "Arial", colour = "black"),
        axis.text     = element_text(size = 6, family = "Arial", colour = "black"),
        axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin   = margin(2, 2, 2, 2)
      )

    # Panel elements only — no title/subtitle/caption, they'd eat most of a
    # tiny figure this size. The stats annotation stays inline (top-left of
    # the panel) instead of a caption, and axis titles are shortened.
    save_figure(p_concordance, file.path(opt$figdir, "04i_cell_line_concordance.pdf"),
                width = 51 / 25.4, height = 31.4 / 25.4)
    message("Figure 4i: Cell-line concordance scatter — done")
  }
} else {
  message("Figures 4b/4h/4i: skipped (Caco2-vs-HEK comparison requires both cell lines; found: ",
          paste(CELL_LINES, collapse = ", "), ")")
}


# ── Figure 3f: Metagene positional profile ────────────────────────────────────
# Where along the transcript do high-confidence sites cluster, aggregated
# across ALL genes of a biotype (not just the top few named ones in 3D) —
# same high-confidence definition as the sequence-context panels (mod_rate
# >= 20%, padj < 0.05), pooled across cell lines since position-along-
# transcript is a property of the RNA, not of the cell line it was measured
# in. Every biotype is shown on its own raw nucleotide position (1..longest
# gene found for that biotype in this dataset) rather than a 0-100% scale —
# miRNA already worked this way (mature-arm position, seed nt 2-8 shaded);
# tRNA/rRNA/snoRNA can too once restricted to their CANONICAL gene_biotype
# (below), which keeps each panel's per-gene length spread tight enough for
# a shared nt axis to be meaningful.
#
# The `biotype` folder label calls are organised under (tRNA/rRNA/snoRNA)
# does not equal the GTF's own `gene_biotype` for every gene_id inside it —
# a handful of gene_ids collide with much longer protein_coding/lncRNA genes
# (max ~900kb), and rRNA further mixes tiny nuclear rRNA/pseudogenes with
# the much longer mitochondrial 12S/16S rRNA. Restricting to the canonical
# GTF biotype for each folder drops those outliers before any nt position
# is computed.
#
# rRNA restricted to true "rRNA" only (dropping "rRNA_pseudogene") 2026-09-
# 07: pseudogenes are genomic copies never incorporated into a functional
# ribosome, so there's no modification literature (e.g. Taoka et al. 2018's
# human rRNA map) to validate or number their positions against -- every
# rRNA site previously surfaced in Figure 4A's top-10 turned out to be a
# pseudogene call. True rRNA (5_8S_rRNA, RNA5S1, RNA5-8SN1, ...) has ample
# qualifying candidates on its own (checked: 3,448 candidate rows clearing
# the usual mod_rate>=0.20/pvalue<0.05/coverage>=50 bar) and, being a fixed,
# non-length-variable transcript per subunit (unlike tRNA), a simple linear
# count from the mature transcript's own 5' end already matches literature
# convention -- no separate universal-numbering system needed here.
CANONICAL_GENE_BIOTYPE <- list(
  tRNA   = "Mt_tRNA",
  rRNA   = "rRNA",
  snoRNA = "snoRNA"
)
is_canonical_biotype <- function(biotype, gene_biotype) {
  dplyr::case_when(
    biotype == "miRNA"  ~ TRUE,
    biotype == "tRNA"   ~ gene_biotype == "Mt_tRNA",
    biotype == "rRNA"   ~ gene_biotype == "rRNA",
    biotype == "snoRNA" ~ gene_biotype == "snoRNA",
    TRUE ~ FALSE
  )
}

# Position normalized to 0-100% of each site's own mature transcript (arm-
# aware for miRNA, gene-span-based for the rest) — the shared axis a
# ridgeline needs to put all four biotypes in ONE panel. Raw-nt axes (used
# in earlier iterations of this figure) can't be shared across biotypes
# whose lengths differ 30x (miRNA ~20-27nt vs snoRNA up to ~791nt).
metagene_sites <- calls %>%
  dplyr::filter(condition == TREAT_COND, mod_rate >= 0.20, padj < 0.05,
                biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA"),
                is_canonical_biotype(biotype, gene_biotype)) %>%
  dplyr::mutate(
    transcript_len = gene_end - gene_start + 1L,
    norm_pos = dplyr::case_when(
      biotype == "miRNA" & !is.na(arm_pos) & arm_len > 1 ~
        (arm_pos - 1) / (arm_len - 1) * 100,
      biotype != "miRNA" & !is.na(mature_pos) & transcript_len > 1 ~
        (mature_pos - 1) / (transcript_len - 1) * 100,
      TRUE ~ NA_real_
    ),
    biotype = factor(biotype, levels = c("snoRNA", "rRNA", "tRNA", "miRNA"))
  ) %>%
  dplyr::filter(!is.na(norm_pos), norm_pos >= 0, norm_pos <= 100)

# Pure theme object (no dependency on metagene_sites) — defined unconditionally
# since Figures 3g/3f-heatmap/etc. below reuse it even when the 3F block itself
# is skipped (e.g. no GTF/miRBase available, so metagene_sites is empty).
metagene_theme <- theme_srnataps() +
  theme(
    text          = element_text(size = 7, family = "Arial", colour = "black"),
    axis.text     = element_text(size = 6, family = "Arial", colour = "black"),
    axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
    plot.title    = element_text(size = 7, family = "Arial", colour = "black"),
    plot.subtitle = element_text(size = 6, family = "Arial", colour = "grey50"),
    panel.grid.minor = element_blank(),
    plot.margin   = margin(4, 8, 4, 4)
  )

if (nrow(metagene_sites) > 0) {
  n_by_bt <- metagene_sites %>% dplyr::count(biotype)
  bt_labels <- setNames(
    paste0(n_by_bt$biotype, " (n=", scales::comma(n_by_bt$n), ")"),
    as.character(n_by_bt$biotype)
  )

  # Back to four SEPARATE density panels, each with its own real numeric
  # y-axis (facet_wrap(scales="free_y")) — the overlaid single-panel version
  # gave a genuinely shared/comparable Density axis, but at the cost of
  # rRNA's high peaks dwarfing miRNA/tRNA's on the same scale. Stacking them
  # instead (facet_wrap, ncol=1) trades cross-biotype comparability for
  # per-biotype shape legibility, same tradeoff as the original per-biotype
  # panels earlier in this session, just now on the shared 0-100% axis.
  metagene_stack <- metagene_sites %>%
    dplyr::mutate(biotype = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")))
  bt_levels_stack <- levels(metagene_stack$biotype)
  n_by_bt_stack <- metagene_stack %>% dplyr::count(biotype)
  bt_labels_stack <- setNames(
    paste0(n_by_bt_stack$biotype, " (n=", scales::comma(n_by_bt_stack$n), ")"),
    as.character(n_by_bt_stack$biotype)
  )

  # Seed (miRNA nt 2-8) and anticodon loop (tRNA nt 32-38) landmarks,
  # converted from literal nt to this panel's normalized 0-100% axis using
  # each biotype's own MEDIAN gene length (every site's norm_pos already
  # accounts for its own gene's length individually; the mark itself needs
  # one representative length to draw a single fixed x-range). Tagged with
  # `biotype` so facet_wrap only renders each band in its own matching
  # facet, not all four.
  mirna_median_len <- stats::median(metagene_stack$arm_len[metagene_stack$biotype == "miRNA"], na.rm = TRUE)
  trna_median_len   <- stats::median(metagene_stack$transcript_len[metagene_stack$biotype == "tRNA"], na.rm = TRUE)
  nt_to_pct <- function(nt, len) (nt - 1) / (len - 1) * 100
  landmarks <- dplyr::tibble(
    biotype = factor(c("miRNA", "tRNA"), levels = bt_levels_stack),
    xmin    = c(nt_to_pct(2, mirna_median_len), nt_to_pct(32, trna_median_len)),
    xmax    = c(nt_to_pct(8, mirna_median_len), nt_to_pct(38, trna_median_len)),
    label   = c("seed", "anticodon")
  )

  p_metagene <- ggplot(metagene_stack, aes(x = norm_pos, colour = biotype, fill = biotype)) +
    geom_rect(
      data = landmarks, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey85", alpha = 0.6, colour = NA
    ) +
    geom_density(adjust = 0.8, alpha = 0.7, linewidth = 0.4) +
    geom_text(
      data = landmarks, aes(x = (xmin + xmax) / 2, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.3, size = 2.11, family = "Arial", colour = "grey30"
    ) +
    annotate("text", x = 1, y = Inf, label = "5'", vjust = 1.4, hjust = -0.1,
             size = 2.46, fontface = "bold", colour = "grey30", family = "Arial") +
    annotate("text", x = 99, y = Inf, label = "3'", vjust = 1.4, hjust = 1.1,
             size = 2.46, fontface = "bold", colour = "grey30", family = "Arial") +
    facet_wrap(~biotype, ncol = 1, scales = "free_y",
              labeller = as_labeller(bt_labels_stack)) +
    scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_fill_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_x_continuous(
      breaks = seq(0, 100, 10), limits = c(0, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title    = "Metagene positional profile of high-confidence sites",
      subtitle = stringr::str_wrap(
        paste0("mod_rate ≥ 20%, padj < 0.05 (", CONDITION_LABELS[TREAT_COND],
              "), pooled Caco2+HEK | position normalized to 0-100% of mature transcript | ",
              "bands = miRNA seed (nt 2-8) & tRNA anticodon loop (nt 32-38), at median gene length | ",
              "y-axis is per-panel (not comparable across biotypes)"),
        width = 100
      ),
      x = "Mature position (%)", y = "Density"
    ) +
    metagene_theme +
    theme(
      strip.background = element_rect(fill = "#1C4062", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = 7, family = "Arial")
    )

  save_figure(p_metagene, file.path(opt$figdir, "04f_metagene_position.pdf"),
              width = 6.3, height = 7.2)
  message("Figure 4f: Metagene positional profile (stacked density panels) — done")
}


# ── Figure 3g: C opportunity vs modification — is the 3F shape just C density? ─
# 3F shows WHERE high-confidence sites cluster (5'→3', normalized position).
# On its own, a peak there is ambiguous: does it mean the assay preferentially
# modifies that region, or are there simply more testable C's there? This
# panel answers that by overlaying, per biotype, the C-OPPORTUNITY density
# (every callable C in the treated condition, canonical-biotype-restricted,
# regardless of whether it was significant) against the SAME modified-site
# density from 3F — both are true densities (integrate to 1), so they share
# one real, comparable y-axis with no scaling trick needed. Where the two
# curves diverge is where positional enrichment is happening independent of
# C availability; where they track together, the "hotspot" is just C density.
c_opportunity_sites <- calls %>%
  dplyr::filter(condition == TREAT_COND,
                biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA"),
                is_canonical_biotype(biotype, gene_biotype)) %>%
  dplyr::mutate(
    transcript_len = gene_end - gene_start + 1L,
    norm_pos = dplyr::case_when(
      biotype == "miRNA" & !is.na(arm_pos) & arm_len > 1 ~
        (arm_pos - 1) / (arm_len - 1) * 100,
      biotype != "miRNA" & !is.na(mature_pos) & transcript_len > 1 ~
        (mature_pos - 1) / (transcript_len - 1) * 100,
      TRUE ~ NA_real_
    ),
    biotype = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA"))
  ) %>%
  dplyr::filter(!is.na(norm_pos), norm_pos >= 0, norm_pos <= 100)

if (nrow(metagene_sites) > 0 && nrow(c_opportunity_sites) > 0) {
  c_context_sites <- metagene_sites %>%
    dplyr::mutate(biotype = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")))

  n_opportunity <- c_opportunity_sites %>% dplyr::count(biotype)
  n_modified    <- c_context_sites %>% dplyr::count(biotype)
  c_context_labels <- setNames(
    paste0(n_opportunity$biotype, "\n", scales::comma(n_modified$n[match(n_opportunity$biotype, n_modified$biotype)]),
          " modified / ", scales::comma(n_opportunity$n), " C's tested"),
    as.character(n_opportunity$biotype)
  )

  mirna_median_len <- stats::median(c_context_sites$arm_len[c_context_sites$biotype == "miRNA"], na.rm = TRUE)
  trna_median_len   <- stats::median(c_context_sites$transcript_len[c_context_sites$biotype == "tRNA"], na.rm = TRUE)
  nt_to_pct <- function(nt, len) (nt - 1) / (len - 1) * 100
  landmarks <- dplyr::tibble(
    biotype = factor(c("miRNA", "tRNA"), levels = c("miRNA", "tRNA", "rRNA", "snoRNA")),
    xmin    = c(nt_to_pct(2, mirna_median_len), nt_to_pct(32, trna_median_len)),
    xmax    = c(nt_to_pct(8, mirna_median_len), nt_to_pct(38, trna_median_len)),
    label   = c("seed", "anticodon")
  )

  # Two-sample Kolmogorov-Smirnov test per biotype: is the modified-site
  # positional distribution actually different from the C-opportunity
  # distribution, or could the two curves plausibly be the same distribution
  # (i.e. the visual "hotspot" is within what C density alone would
  # produce)? D = max distance between the two empirical CDFs (0 = identical
  # distributions, larger = more different); p from the asymptotic KS null.
  ks_stats <- dplyr::bind_rows(lapply(levels(c_context_sites$biotype), function(bt) {
    kt <- suppressWarnings(stats::ks.test(
      c_context_sites$norm_pos[c_context_sites$biotype == bt],
      c_opportunity_sites$norm_pos[c_opportunity_sites$biotype == bt]
    ))
    dplyr::tibble(
      biotype = factor(bt, levels = levels(c_context_sites$biotype)),
      label = paste0("KS D=", sprintf("%.2f", kt$statistic),
                    ", p", if (kt$p.value < 0.001) "<0.001" else paste0("=", sprintf("%.3f", kt$p.value)))
    )
  }))

  # Both series are true densities (integrate to 1), so they share one real,
  # comparable y-axis with no scaling trick needed.
  p_c_context <- ggplot() +
    geom_rect(
      data = landmarks, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey85", alpha = 0.6, colour = NA
    ) +
    geom_density(
      data = c_opportunity_sites, aes(x = norm_pos),
      fill = "grey55", colour = "grey40", alpha = 0.35, linewidth = 0.35, adjust = 0.8
    ) +
    geom_density(
      data = c_context_sites, aes(x = norm_pos, fill = biotype, colour = biotype),
      alpha = 0.55, linewidth = 0.4, adjust = 0.8
    ) +
    geom_text(
      data = landmarks, aes(x = (xmin + xmax) / 2, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.3, size = 2.11, family = "Arial", colour = "grey30"
    ) +
    geom_text(
      data = ks_stats, aes(x = 98, y = Inf, label = label),
      inherit.aes = FALSE, hjust = 1, vjust = 1.3, size = 2.11, family = "Arial",
      colour = "grey30", fontface = "italic"
    ) +
    facet_wrap(~biotype, ncol = 1, scales = "free_y",
              labeller = as_labeller(c_context_labels)) +
    scale_fill_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_x_continuous(
      breaks = seq(0, 100, 10), limits = c(0, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title    = "C opportunity vs modification: is the positional pattern just C density?",
      subtitle = stringr::str_wrap(
        paste0("Grey = every callable C, ", CONDITION_LABELS[TREAT_COND], " (C opportunity) | colour = high-confidence ",
              "modified sites (same as Fig. 3F) | mod_rate ≥ 20%, padj < 0.05 | pooled Caco2+HEK | ",
              "KS D/p = two-sample Kolmogorov-Smirnov test, modified vs opportunity distribution (D=0 identical, ",
              "larger D/small p = positional pattern is NOT just explained by C density) | ",
              "y-axis is per-panel (not comparable across biotypes)"),
        width = 100
      ),
      x = "Mature position (%)", y = "Density"
    ) +
    metagene_theme +
    theme(
      strip.background = element_rect(fill = "#1C4062", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = 6, family = "Arial", lineheight = 0.95)
    )

  save_figure(p_c_context, file.path(opt$figdir, "04g_c_opportunity_context.pdf"),
              width = 6.3, height = 7.2)
  message("Figure 4g: C opportunity vs modification context — done")

  # ── Figure 3j: CpG vs non-CpG — where does the reported m5C signal sit? ────
  # Same two universes as 3G (C opportunity = every callable C; modified =
  # high-confidence sites), now split by dinucleotide context instead of
  # position: is the downstream (3') base a G (CpG) or not (non-CpG)? Paired
  # bars per biotype so a CpG enrichment in the modified set can be checked
  # against the CpG content of the underlying testable-C universe, not just
  # eyeballed against a "50/50 would be unbiased" assumption — the C
  # opportunity universe itself is very unlikely to be 50/50 CpG.
  classify_cpg <- function(df) {
    df %>%
      dplyr::filter(nchar(context) == 3) %>%
      dplyr::mutate(cpg = dplyr::if_else(toupper(substr(context, 3, 3)) == "G", "CpG", "non-CpG"))
  }
  cpg_opportunity <- classify_cpg(c_opportunity_sites) %>%
    dplyr::count(biotype, cpg) %>%
    dplyr::group_by(biotype) %>%
    dplyr::mutate(pct = n / sum(n), universe = "C opportunity") %>%
    dplyr::ungroup()
  cpg_modified <- classify_cpg(c_context_sites) %>%
    dplyr::count(biotype, cpg) %>%
    dplyr::group_by(biotype) %>%
    dplyr::mutate(pct = n / sum(n), universe = "Modified (HC)") %>%
    dplyr::ungroup()
  cpg_all <- dplyr::bind_rows(cpg_opportunity, cpg_modified) %>%
    dplyr::mutate(
      biotype  = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")),
      universe = factor(universe, levels = c("C opportunity", "Modified (HC)")),
      cpg      = factor(cpg, levels = c("non-CpG", "CpG"))
    )

  # Fisher's exact test per biotype: is the modified set's CpG fraction
  # different from the opportunity universe's CpG fraction, beyond what's
  # expected by chance given each biotype's own sample sizes?
  cpg_tests <- dplyr::bind_rows(lapply(levels(cpg_all$biotype), function(bt) {
    opp <- cpg_opportunity %>% dplyr::filter(biotype == bt)
    mod <- cpg_modified %>% dplyr::filter(biotype == bt)
    tbl <- matrix(
      c(mod$n[mod$cpg == "CpG"], mod$n[mod$cpg == "non-CpG"],
        opp$n[opp$cpg == "CpG"], opp$n[opp$cpg == "non-CpG"]),
      nrow = 2, dimnames = list(c("CpG", "non-CpG"), c("Modified", "Opportunity"))
    )
    ft <- stats::fisher.test(tbl)
    dplyr::tibble(
      biotype = factor(bt, levels = levels(cpg_all$biotype)),
      label = paste0("OR=", sprintf("%.2f", ft$estimate),
                    "\np", if (ft$p.value < 0.001) "<0.001" else paste0("=", sprintf("%.3f", ft$p.value)))
    )
  }))

  CPG_COLOURS <- c("CpG" = "#2E7D32", "non-CpG" = "#B0BEC5")

  p_cpg <- ggplot(cpg_all, aes(x = universe, y = pct, fill = cpg)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
    geom_text(
      data = dplyr::filter(cpg_all, cpg == "CpG"),
      aes(label = scales::percent(pct, accuracy = 1)),
      position = position_stack(vjust = 0.5), size = 2.11, family = "Arial", colour = "white", fontface = "bold"
    ) +
    geom_text(
      data = cpg_tests, aes(x = 1.5, y = 1.09, label = label),
      inherit.aes = FALSE, size = 1.76, family = "Arial", colour = "grey30", fontface = "italic", lineheight = 0.85
    ) +
    facet_wrap(~biotype, nrow = 1) +
    scale_fill_manual(values = CPG_COLOURS, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.22),
                       breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "CpG vs non-CpG context of reported m5C sites",
      subtitle = stringr::str_wrap(
        paste0("% of sites where the downstream (3') base is G (CpG) vs not (non-CpG), C opportunity (every callable C) ",
              "vs high-confidence modified sites (same criteria as Fig. 3G) | mod_rate ≥ 20%, padj < 0.05 | ",
              CONDITION_LABELS[TREAT_COND], ", pooled Caco2+HEK | OR/p = Fisher's exact test, modified vs opportunity CpG odds"),
        width = 58
      ),
      x = NULL, y = "% of sites"
    ) +
    metagene_theme +
    theme(
      strip.background = element_rect(fill = "#1C4062", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = 6, family = "Arial"),
      axis.text.x       = element_text(size = 6, family = "Arial", colour = "black", angle = 20, hjust = 1),
      legend.text       = element_text(size = 7, family = "Arial", colour = "black"),
      legend.position   = "bottom"
    )

  # Fix each facet panel box to an exact physical size (18.8 x 13.1mm),
  # independent of strip/axis/legend sizing, by editing the built gtable
  # directly — ggsave's width/height only control the overall device.
  gt_cpg     <- ggplotGrob(p_cpg)
  is_panel   <- grepl("^panel($|-)", gt_cpg$layout$name)
  panel_row  <- unique(gt_cpg$layout$t[is_panel])
  panel_col  <- unique(gt_cpg$layout$l[is_panel])
  gt_cpg$widths[panel_col]  <- unit(18.8, "mm")
  gt_cpg$heights[panel_row] <- unit(13.1, "mm")

  pdf(NULL)
  total_w_in <- sum(grid::convertWidth(gt_cpg$widths, "in", valueOnly = TRUE))
  total_h_in <- sum(grid::convertHeight(gt_cpg$heights, "in", valueOnly = TRUE))
  dev.off()

  base_out_j <- file.path(opt$figdir, "04j_cpg_context")
  ggsave(paste0(base_out_j, ".pdf"), gt_cpg, width = total_w_in, height = total_h_in,
        device = cairo_pdf, family = "Arial")
  ggsave(paste0(base_out_j, ".png"), gt_cpg, width = total_w_in, height = total_h_in, dpi = 300)
  ggsave(paste0(base_out_j, ".svg"), gt_cpg, width = total_w_in, height = total_h_in, device = "svg")
  message("  Saved: ", paste0(base_out_j, ".pdf"), " | each panel fixed at 18.8x13.1mm")
  message("Figure 4j: CpG vs non-CpG context — done")
}

# ── Figure 3k: Top 10 modified sites per biotype (Caco2 vs HEK) ─────────────
# Same high-confidence criteria as every other site-level Fig. 3 panel
# (mod_rate >= 20%, pvalue < 0.05, canonical biotype only, TET+PB treatment)
# -- ranks each biotype's own sites by median modification rate (across cell
# lines) and shows the top 10, split by cell line (shape) with point size
# encoding coverage depth.
#
# TOP10_MIN_COVERAGE guards against a small-N artifact discovered 2026-09-06:
# with only the pipeline's usual min-cov=5 floor, "top by raw rate" is won
# almost entirely by sites with 5-20 total reads that happen to be 100%
# modified by chance -- e.g. every one of the (pre-fix) top-10 miRNA/rRNA/
# snoRNA sites had coverage <=20 and rate exactly 1.0, while tRNA's top-10
# (naturally >=35x deeper, mitochondrial) looked biologically sane. 50 reads
# clears that failure mode while still leaving hundreds of eligible sites per
# biotype (checked: >=89 for the shallowest biotype, tRNA) -- comfortably
# enough for a top-10 (or top-5, rRNA) selection.
TOP10_MIN_COVERAGE <- 50
hc_calls_top10 <- calls %>%
  dplyr::filter(condition == TREAT_COND, mod_rate >= 0.20, pvalue < 0.05,
                coverage >= TOP10_MIN_COVERAGE,
                biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA"),
                is_canonical_biotype(biotype, gene_biotype)) %>%
  dplyr::mutate(
    transcript_len = gene_end - gene_start + 1L,
    raw_pos  = dplyr::if_else(biotype == "miRNA", arm_pos, mature_pos),
    norm_pos = dplyr::case_when(
      biotype == "miRNA" & !is.na(arm_pos) & arm_len > 1 ~ (arm_pos - 1) / (arm_len - 1) * 100,
      biotype != "miRNA" & !is.na(mature_pos) & transcript_len > 1 ~ (mature_pos - 1) / (transcript_len - 1) * 100,
      TRUE ~ NA_real_
    ),
    # tRNA's own numeric part of the Sprinzl label (e.g. "61A" -> 61), used
    # only for the anticodon-loop landmark range check below -- a linear
    # mature_pos was never the right thing to compare against 32-38 (see
    # sprinzl_position's own definition above), and this restores that
    # literature-defined landmark to the numbering system it's actually
    # defined in.
    sprinzl_numeric = dplyr::if_else(
      biotype == "tRNA" & !is.na(sprinzl_position),
      as.numeric(gsub("[^0-9]", "", sprinzl_position)),
      NA_real_
    ),
    # Same landmark definitions as Fig. 3f/3g: miRNA seed (arm nt2-8), tRNA
    # anticodon loop (Sprinzl 32-38). No defined landmark for rRNA/snoRNA.
    landmark = dplyr::case_when(
      biotype == "miRNA" & !is.na(arm_pos) & arm_pos >= 2 & arm_pos <= 8 ~ "seed",
      biotype == "tRNA" & !is.na(sprinzl_numeric) & sprinzl_numeric >= 32 & sprinzl_numeric <= 38 ~ "anticodon",
      TRUE ~ NA_character_
    ),
    # Display label: tRNA shows Sprinzl/universal numbering (a different
    # system than the x-axis position below, which stays linear -- Sprinzl
    # slots aren't evenly spaced, so they can't drive point placement, only
    # the printed label). miRNA/rRNA/snoRNA keep their existing linear "ntN"
    # label, tagged separately from tRNA's "UN" so the two numbering
    # systems are never visually ambiguous in the same panel.
    #
    # miRNA precursor fallback: a site can be a real, high-confidence call
    # that sits outside every annotated mature 5p/3p arm interval (e.g. in
    # the few nucleotides of hairpin between the precursor's own 5' end and
    # where Drosha/Dicer processing actually begins the mature product) --
    # arm_pos is then NA not because of a data problem but because there
    # genuinely is no mature-arm position for that site. Rather than leave
    # the label blank, fall back to the PRECURSOR-relative linear position
    # (`mature_pos`, already computed for every biotype above) with its own
    # "pre" prefix -- explicitly marked as pri-miRNA/precursor coordinates,
    # never mistaken for a mature-arm position.
    label_kind = dplyr::case_when(
      biotype == "tRNA" ~ "sprinzl",
      biotype == "miRNA" & !is.na(arm_pos) ~ "mature_arm",
      biotype == "miRNA" & is.na(arm_pos) & !is.na(mature_pos) & mature_pos >= 1 & transcript_len > 1 ~ "precursor",
      biotype %in% c("rRNA", "snoRNA") & !is.na(mature_pos) ~ "mature_linear",
      TRUE ~ NA_character_
    ),
    label_pos = dplyr::case_when(
      label_kind == "sprinzl"       ~ sprinzl_position,
      label_kind == "mature_arm"    ~ as.character(arm_pos),
      label_kind == "precursor"     ~ as.character(mature_pos),
      label_kind == "mature_linear" ~ as.character(mature_pos),
      TRUE ~ NA_character_
    )
  ) %>%
  # Gene-span boundary edge cases (site start falling just outside the
  # annotated [gene_start, gene_end], same coordinate-system artifact
  # already guarded against for the metagene panels via
  # `filter(norm_pos >= 0, norm_pos <= 100)`) can push a handful of sites to
  # nt0/negative % or beyond 100% -- blank the position display for those
  # rather than showing an impossible value; the site itself still counts
  # for ranking, only its position label is suppressed. This guards
  # raw_pos/norm_pos/landmark only -- label_pos has its own validity baked
  # into label_kind above (deliberately independent, since a miRNA
  # "precursor" label is expected to have norm_pos/landmark = NA; that's not
  # an error state for the label itself).
  dplyr::mutate(
    valid_pos = !is.na(norm_pos) & norm_pos >= 0 & norm_pos <= 100 & !is.na(raw_pos) & raw_pos >= 1,
    raw_pos   = dplyr::if_else(valid_pos, as.numeric(raw_pos), NA_real_),
    norm_pos  = dplyr::if_else(valid_pos, norm_pos, NA_real_),
    landmark  = dplyr::if_else(valid_pos, landmark, NA_character_)
  )

if (nrow(hc_calls_top10) > 0) {
  site_cell_top10 <- hc_calls_top10 %>%
    dplyr::mutate(site_id = paste0(biotype, ":", chrom, ":", start)) %>%
    dplyr::group_by(biotype, site_id, site_label, cell_line) %>%
    dplyr::summarise(
      mod_rate = mean(mod_rate), coverage = mean(coverage),
      raw_pos = dplyr::first(raw_pos), norm_pos = dplyr::first(norm_pos), landmark = dplyr::first(landmark),
      label_pos = dplyr::first(label_pos), label_kind = dplyr::first(label_kind),
      .groups = "drop"
    )

  # A gene can carry more than one high-confidence C position (e.g. MT-TS2,
  # SNORD43) -- site_id is position-level but site_label is gene-level, so
  # multiple top-10 rows can share one display label. Append the mature
  # position (nt / % along the transcript, same norm_pos metric as Fig.
  # 3f/3g) and landmark flag to every label -- this both surfaces the
  # positional context the user wants AND naturally disambiguates same-gene
  # positions in almost all cases. A genomic-coordinate suffix remains as a
  # last-resort tie-breaker for the rare case mature position is still NA
  # or identical (e.g. unmapped position, or a duplicate row).
  site_rank_top10 <- site_cell_top10 %>%
    dplyr::group_by(biotype, site_id, site_label) %>%
    dplyr::summarise(
      median_rate = stats::median(mod_rate),
      raw_pos = dplyr::first(raw_pos), norm_pos = dplyr::first(norm_pos), landmark = dplyr::first(landmark),
      label_pos = dplyr::first(label_pos), label_kind = dplyr::first(label_kind),
      .groups = "drop"
    ) %>%
    dplyr::group_by(biotype) %>%
    # rRNA gets top 5 (not 10) -- its rates cluster so tightly (mostly
    # 85-100%) that a full 10 crowds the panel with near-duplicate points;
    # every other biotype keeps the usual top 10. slice_max's n must be a
    # constant, so this needs two separate calls rather than a per-group
    # expression.
    dplyr::group_split() %>%
    lapply(function(df) {
      n_top <- if (df$biotype[1] == "rRNA") 5L else 10L
      dplyr::slice_max(df, median_rate, n = n_top, with_ties = FALSE)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(biotype, site_label) %>%
    dplyr::mutate(
      site_label = dplyr::if_else(
        rep(dplyr::n() > 1, dplyr::n()),
        paste0(site_label, " [", sub("^[^:]+:", "", site_id), "]"),
        site_label
      )
    ) %>%
    dplyr::ungroup()

  # Row order: a gene with >1 top-10 site (e.g. MT-TS2, SNORD43) previously
  # got its rows scattered across the panel, interleaved with other genes,
  # because ordering was purely by each individual site's own rate. Group in
  # two nested levels: (1) paralog FAMILY (e.g. SNORD62A and SNORD62B both
  # -> "SNORD62"), so related genes land next to each other too, not just
  # repeat sites of the exact same gene; (2) the exact gene within that
  # family. Both levels are ordered by their own best rate (so the overall
  # "highest near the top" shape is kept), with each site's own rate as the
  # final tie-breaker -- every row for a family/gene ends up adjacent.
  #
  # The family regex strips a trailing letter(+digit) ONLY when it follows a
  # digit (lookbehind) -- this catches true paralog suffixes (SNORD62A/B,
  # MIR200A/B/C) without misfiring on tRNA gene symbols like "MT-TM", where
  # the trailing letter is the amino-acid code, not a paralog suffix, and
  # isn't preceded by a digit.
  gene_order_top10 <- site_rank_top10 %>%
    dplyr::mutate(
      gene_base   = sub(" \\[.*\\]$", "", site_label),
      family_base = sub("(?<=[0-9])[A-Z][0-9]*$", "", gene_base, perl = TRUE)
    ) %>%
    dplyr::group_by(biotype, family_base) %>%
    dplyr::mutate(family_best_rate = max(median_rate)) %>%
    dplyr::group_by(biotype, gene_base) %>%
    dplyr::mutate(gene_best_rate = max(median_rate)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(biotype, family_best_rate, family_base, gene_best_rate, gene_base, median_rate) %>%
    dplyr::pull(site_label) %>%
    unique()

  top10_sites <- site_cell_top10 %>%
    dplyr::select(-site_label) %>%
    dplyr::inner_join(
      dplyr::select(site_rank_top10, biotype, site_id, site_label, median_rate),
      by = c("biotype", "site_id")
    ) %>%
    dplyr::mutate(
      biotype    = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")),
      cell_line  = factor(cell_line, levels = CELL_LINES),
      site_label = factor(site_label, levels = gene_order_top10),
      # Position label, shown directly above each point instead of baked
      # into the y-axis text. tRNA uses Sprinzl/universal tRNA numbering
      # ("U" prefix), a different system from every other biotype's linear
      # mature-transcript position ("nt" prefix) -- prefixed distinctly so
      # the two are never mistaken for the same numbering (see label_pos's
      # own definition above for why that mistake is easy to make otherwise).
      # "U" = Sprinzl/universal tRNA position; "nt" = linear mature position
      # (mature arm for miRNA, mature transcript for rRNA/snoRNA); "pre" =
      # miRNA precursor/pri-miRNA position -- used only when a site falls
      # outside every annotated mature arm, so there's no mature position to
      # report at all (see label_kind's own definition above).
      pos_text = dplyr::case_when(
        is.na(label_pos) ~ NA_character_,
        label_kind == "sprinzl"   ~ paste0("U", label_pos),
        label_kind == "precursor" ~ paste0("pre", label_pos),
        TRUE ~ paste0("nt", label_pos)
      ),
      # rRNA's rates cluster so tightly that Caco2/HEK points for the same
      # site often sit almost exactly on top of each other horizontally --
      # nudge them vertically apart (stacked within the row) there only;
      # every other biotype keeps both points on the row's centre line.
      y_nudge = dplyr::if_else(biotype == "rRNA",
                               dplyr::if_else(cell_line == "Caco2", 0.18, -0.18), 0)
    )

  # Position is a property of the SITE, not of the cell-line point -- it's
  # identical for a site's Caco2 and HEK rows. Labelling both would just
  # duplicate the same text twice at nearly the same spot, so anchor one
  # label per site above whichever point has the higher rate (usually the
  # one with more clear space to its right).
  # Anchor at the point with the LOWER rate of the pair (leftmost) so the
  # label -- growing further left, hjust=1 -- has room to sit clear of both
  # cell-line points instead of colliding with whichever one is further
  # right. Flips to grow rightward instead (hjust=0) when that point is
  # already close to the x=40 panel edge, so the label doesn't clip off the
  # left side (happened to MT-TF [MT:645]'s "nt69" at x~48).
  pos_labels_top10 <- top10_sites %>%
    dplyr::filter(!is.na(pos_text)) %>%
    dplyr::group_by(biotype, site_label) %>%
    dplyr::slice_min(mod_rate, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      near_left_edge = mod_rate * 100 < 50,
      label_x     = dplyr::if_else(near_left_edge, mod_rate * 100 + 2, mod_rate * 100 - 2),
      label_hjust = dplyr::if_else(near_left_edge, 0, 1)
    )

  p_top10 <- ggplot(top10_sites, aes(x = mod_rate * 100, y = site_label)) +
    geom_point(aes(colour = biotype, shape = cell_line, size = log10(coverage)),
              alpha = 0.9, stroke = 0.6,
              position = ggplot2::position_nudge(y = top10_sites$y_nudge)) +
    geom_text(
      data = pos_labels_top10, aes(label = pos_text, x = label_x, hjust = label_hjust),
      vjust = 0.35, size = 1.76,
      family = "Arial", colour = "grey30", na.rm = TRUE, show.legend = FALSE,
      position = ggplot2::position_nudge(y = pos_labels_top10$y_nudge)
    ) +
    facet_wrap(~biotype, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_shape_manual(values = CELL_SHAPES, name = "Cell line") +
    # Full 1-4 legend key range requested 2026-09-07 (previously trimmed to
    # 2 keys only for space) -- wrapped to 2 rows via the size guide so the
    # 4 keys don't clip off the narrow single-column device's right edge.
    # `limits` explicitly forced to c(1, 4): the actual data never goes
    # below log10(coverage)~1.7 (TOP10_MIN_COVERAGE=50), so break=1 would
    # otherwise silently disappear from the legend (ggplot drops breaks
    # outside the scale's data range) even though the request is for a
    # legend spanning the full 1-4 axis, not just the range this panel's
    # data happens to occupy.
    scale_size_continuous(name = "log10(coverage)", range = c(1, 4), breaks = c(1, 2, 3, 4),
                          limits = c(1, 4), guide = guide_legend(nrow = 2)) +
    scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 20),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      # Wrapped: the single-column stack is only ~35mm+axis wide now, far
      # narrower than the old 2-column layout these were sized for -- an
      # unwrapped title/subtitle just ran off the edge of the device.
      title    = stringr::str_wrap(
        paste0("Top 10 modified sites per biotype (", CONDITION_LABELS[TREAT_COND], ")"), width = 30),
      subtitle = stringr::str_wrap(
        paste0(paste(CELL_LINES, collapse = " and "), " — ranked by median modification rate"), width = 30),
      # tRNA positions are a genuinely different numbering system from every
      # other biotype's (Sprinzl/universal tRNA numbering, not a linear
      # count from the gene's own start) -- spelled out here since the two
      # look identical at a glance ("U49" vs "nt49") but are not comparable
      # numbers at all.
      caption  = stringr::str_wrap(
        paste0("U = Sprinzl/universal tRNA position; nt = linear mature-transcript position ",
               "(mature arm for miRNA, mature transcript for rRNA/snoRNA); pre = miRNA precursor ",
               "position (site falls outside every annotated mature arm)"),
        width = 34),
      x = "Modification rate (%)", y = NULL
    ) +
    metagene_theme +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text        = element_text(colour = "black", face = "bold", size = 6, family = "Arial"),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.text       = element_text(size = 7, family = "Arial", colour = "black"),
      legend.title      = element_text(size = 7, family = "Arial", colour = "black"),
      plot.caption      = element_text(size = 5.5, family = "Arial", colour = "grey50", hjust = 0)
    )

  # 4 panels stacked in one column (miRNA/tRNA/rRNA/snoRNA top to bottom),
  # each fixed to exactly 35mm wide, independent of axis/strip/legend
  # sizing -- ggsave's width/height only control the overall device. Point
  # size and every font size are untouched by this; only the gtable's
  # panel-cell width/height units change. Height is 45mm for the 10-row
  # panels, but rRNA only has 5 rows (top 5, not top 10) -- scaling its
  # height to 22.5mm keeps the same per-row spacing as the others instead
  # of leaving it with tall, sparse blank space.
  gt_top10   <- ggplotGrob(p_top10)
  is_panel_k <- grepl("^panel($|-)", gt_top10$layout$name)
  panel_row_k <- unique(gt_top10$layout$t[is_panel_k])
  panel_col_k <- unique(gt_top10$layout$l[is_panel_k])
  gt_top10$widths[panel_col_k]  <- unit(35, "mm")
  # panel_row_k is sorted top-to-bottom, matching the biotype factor order
  # (miRNA, tRNA, rRNA, snoRNA) that ncol=1 stacks in.
  gt_top10$heights[panel_row_k] <- unit(c(45, 45, 22.5, 45), "mm")

  pdf(NULL)
  total_w_in_k <- sum(grid::convertWidth(gt_top10$widths, "in", valueOnly = TRUE))
  total_h_in_k <- sum(grid::convertHeight(gt_top10$heights, "in", valueOnly = TRUE))
  dev.off()

  base_out_k <- file.path(opt$figdir, "04k_top10_sites_per_biotype")
  ggsave(paste0(base_out_k, ".pdf"), gt_top10, width = total_w_in_k, height = total_h_in_k,
        device = cairo_pdf, family = "Arial")
  ggsave(paste0(base_out_k, ".png"), gt_top10, width = total_w_in_k, height = total_h_in_k, dpi = 300)
  ggsave(paste0(base_out_k, ".svg"), gt_top10, width = total_w_in_k, height = total_h_in_k, device = "svg")
  message("  Saved: ", paste0(base_out_k, ".pdf"), " | panels fixed at 35mm wide, 45mm tall (22.5mm for rRNA)")
  message("Figure 4k: Top 10 modified sites per biotype — done")

  # Heatmap alternative: site x cell_line grid, fill = modification rate —
  # same white->TREAT_COND gradient convention as Fig. 3E's trinucleotide
  # context heatmap. tidyr::expand_grid ensures a blank (grey) tile is drawn
  # for any site missing one cell line, rather than silently omitting it.
  heat_grid_top10 <- tidyr::expand_grid(
    dplyr::select(site_rank_top10, biotype, site_id, site_label, median_rate),
    cell_line = CELL_LINES
  ) %>%
    dplyr::left_join(
      dplyr::select(site_cell_top10, biotype, site_id, cell_line, mod_rate, coverage),
      by = c("biotype", "site_id", "cell_line")
    ) %>%
    dplyr::mutate(
      biotype    = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")),
      cell_line  = factor(cell_line, levels = CELL_LINES),
      # Same gene-grouped row order as the dot-plot version above.
      site_label = factor(site_label, levels = gene_order_top10),
      text_col   = dplyr::if_else(!is.na(mod_rate) & mod_rate > 0.55, "white", "black")
    )

  p_top10_heat <- ggplot(heat_grid_top10, aes(x = cell_line, y = site_label, fill = mod_rate)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(
      aes(label = dplyr::if_else(is.na(mod_rate), "", paste0(round(mod_rate * 100), "%")), colour = text_col),
      size = 2.11, family = "Arial"
    ) +
    scale_colour_identity() +
    facet_wrap(~biotype, ncol = 2, scales = "free_y") +
    scale_fill_gradient(
      low = "#FFFFFF", high = CONDITION_COLOURS[TREAT_COND], na.value = "grey92",
      limits = c(0, 1), name = "Mod\nrate", labels = function(x) paste0(round(x * 100), "%")
    ) +
    labs(
      title    = paste0("Top 10 modified sites per biotype (", CONDITION_LABELS[TREAT_COND], ")"),
      subtitle = "ranked by median modification rate across cell lines",
      x = NULL, y = NULL
    ) +
    metagene_theme +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text        = element_text(colour = "black", face = "bold", size = 6, family = "Arial"),
      panel.grid        = element_blank(),
      axis.text.x       = element_text(size = 6, family = "Arial", colour = "black"),
      legend.text       = element_text(size = 7, family = "Arial", colour = "black"),
      legend.title      = element_text(size = 7, family = "Arial", colour = "black")
    )

  save_figure(p_top10_heat, file.path(opt$figdir, "04k_top10_sites_heatmap.pdf"), width = 6.5, height = 6.2)
  message("Figure 4k (heatmap alt): Top 10 modified sites per biotype — done")
}


# ── Figure 3c: Condition comparison scatter (site level) ─────────────────────
cond_wide <- calls %>%
  dplyr::filter(condition %in% c(CTRL_COND, TREAT_COND)) %>%
  dplyr::mutate(site_id = paste0(chrom, ":", start)) %>%
  dplyr::group_by(biotype, site_id, site_label, cell_line, condition) %>%
  dplyr::summarise(mod_rate = median(mod_rate), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = mod_rate) %>%
  dplyr::filter(!is.na(.data[[CTRL_COND]]), !is.na(.data[[TREAT_COND]]))

if (nrow(cond_wide) > 0) {
  # Rows = cell line (Caco2 top, HEK bottom), columns = biotype — the same
  # facet_grid(cell_line ~ biotype) arrangement as Figure 3A, for visual
  # consistency between the two panels.
  cond_wide <- cond_wide %>%
    dplyr::mutate(
      above_diag = .data[[TREAT_COND]] > .data[[CTRL_COND]],
      biotype    = factor(biotype, levels = BIOTYPES_PLOT),
      cell_line  = factor(cell_line, levels = c("Caco2", "HEK"))
    )

  # Colour encodes which condition dominates at a site, using the same
  # fixed condition palette as Figure 3A.
  DOMINANT_COLOURS <- c(
    "PB+ (TET-) only" = COND_PB_CTRL_COLOUR,
    "PB+ & TET+"       = COND_TREAT_COLOUR
  )
  cond_wide <- cond_wide %>%
    dplyr::mutate(dominant = dplyr::if_else(
      above_diag, "PB+ & TET+", "PB+ (TET-) only"
    ))

  # RMSE from the identity line (root-mean-square of TET+PB minus PB-only, in
  # percentage points) captures both the systematic TET-induced shift and the
  # scatter around it in one number per panel; Pearson r is added alongside as
  # a check that the two conditions still track each other, not just diverge.
  scatter_stats <- cond_wide %>%
    dplyr::group_by(biotype, cell_line) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      rmse = sqrt(mean((.data[[TREAT_COND]] - .data[[CTRL_COND]])^2)) * 100,
      r    = stats::cor(.data[[CTRL_COND]], .data[[TREAT_COND]]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(label = sprintf("RMSE = %.1fpp\nr = %.2f", rmse, r))

  p_scatter <- ggplot(cond_wide,
                      aes(x = .data[[CTRL_COND]], y = .data[[TREAT_COND]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.4) +
    geom_point(aes(colour = dominant), shape = 16, alpha = 0.3, size = 0.9) +
    geom_text(data = scatter_stats,
              aes(x = 0.03, y = 0.97, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 1.85,
              lineheight = 0.88, family = "Arial", colour = "#1C4062") +
    facet_grid(cell_line ~ biotype) +
    scale_colour_manual(values = DOMINANT_COLOURS, name = "Dominant condition") +
    scale_x_continuous(labels = function(x) x * 100, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(labels = function(y) y * 100, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(
      title    = "Condition comparison: PB-only vs TET+PB",
      subtitle = "Caco2 (top) and HEK (bottom), per biotype | colour = dominant condition",
      x        = paste0(CONDITION_LABELS[CTRL_COND],  " mod rate (%)"),
      y        = paste0(CONDITION_LABELS[TREAT_COND], " mod rate (%)"),
      caption  = "Dashed line = identity"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    theme_srnataps() +
    theme(
      legend.position = "bottom",
      text          = element_text(size = 7, family = "Arial", colour = "black"),
      axis.text     = element_text(size = 7, family = "Arial", colour = "black"),
      axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
      plot.title    = element_text(size = 7, family = "Arial", colour = "black"),
      plot.subtitle = element_text(size = 6, family = "Arial", colour = "grey50"),
      legend.text   = element_text(size = 7, family = "Arial", colour = "black"),
      legend.title  = element_text(size = 7, family = "Arial", colour = "black"),
      strip.text    = element_text(size = 7, family = "Arial", colour = "white", face = "bold"),
      plot.caption  = element_text(size = 6, family = "Arial", colour = "grey50")
    )
  # Fix each facet panel to an exact physical width (29mm), independent of
  # axis/legend/strip sizing, by editing the built gtable directly — ggsave's
  # width/height only control the overall device, not individual panel size.
  gt_scatter  <- ggplotGrob(p_scatter)
  panel_cols  <- unique(gt_scatter$layout$l[grepl("^panel-", gt_scatter$layout$name)])
  gt_scatter$widths[panel_cols] <- unit(29, "mm")

  # Only the panel columns were fixed to an absolute unit; row heights are
  # still "null" (flexible) and can't be measured outside a real render, so
  # height stays an explicit target (matching Figure 3A) rather than derived.
  pdf(NULL)
  total_w_in <- sum(grid::convertWidth(gt_scatter$widths, "in", valueOnly = TRUE))
  dev.off()
  total_h_in <- 3.6

  base_out <- file.path(opt$figdir, "04c_condition_comparison")
  ggsave(paste0(base_out, ".pdf"), gt_scatter, width = total_w_in, height = total_h_in,
        device = cairo_pdf, family = "Arial")
  ggsave(paste0(base_out, ".png"), gt_scatter, width = total_w_in, height = total_h_in, dpi = 300)
  ggsave(paste0(base_out, ".svg"), gt_scatter, width = total_w_in, height = total_h_in, device = "svg")
  message("  Saved: ", paste0(base_out, ".pdf"), " | panel width fixed at 29mm")
  message("Figure 4c: Condition comparison scatter — done")
}

# ── Figure 3d: Waterfall plot (sites ranked by mod_rate) ─────────────────────
waterfall_data <- calls %>%
  dplyr::filter(condition == TREAT_COND) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::arrange(desc(mod_rate), .by_group = TRUE) %>%
  dplyr::mutate(rank = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cell_line) %>%
  dplyr::slice_head(n = 500) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(biotype = factor(biotype, levels = names(BIOTYPE_COLOURS)))

if (nrow(waterfall_data) > 0) {
  p_waterfall <- ggplot(waterfall_data,
                        aes(x = rank, y = mod_rate * 100,
                            colour = biotype, shape = biotype)) +
    geom_point(size = 1.5, alpha = 0.85) +
    facet_wrap(~ cell_line) +
    scale_colour_manual(values = BIOTYPE_COLOURS, name = "Biotype") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 3, 4, 1), name = "Biotype") +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    scale_x_continuous(breaks = seq(0, 500, by = 100)) +
    labs(
      title    = paste0("Top 500 modified sites (", CONDITION_LABELS[TREAT_COND], ")"),
      subtitle = "Sites ranked by modification rate | coloured by RNA biotype",
      x        = "Rank",
      y        = "Modification rate (%)"
    ) +
    theme_srnataps() +
    theme(legend.position = "bottom")
  save_figure(p_waterfall, file.path(opt$figdir, "04d_waterfall.pdf"),
              width = 9, height = 5)
  message("Figure 4d: Waterfall plot — done")
}

# ── Figure 3e: Trinucleotide context heatmap ─────────────────────────────────

context_data <- calls %>%
  dplyr::filter(condition == TREAT_COND, nchar(context) == 3) %>%
  dplyr::mutate(
    ctx_up   = substr(context, 1, 1),
    ctx_down = substr(context, 3, 3)
  ) %>%
  dplyr::group_by(biotype, cell_line, ctx_up, ctx_down) %>%
  dplyr::summarise(
    mean_mod = mean(mod_rate),
    n_sites  = n(),
    .groups  = "drop"
  ) %>%
  dplyr::filter(n_sites >= 5)

if (nrow(context_data) > 0) {
  context_data <- context_data %>%
    dplyr::mutate(text_col = ifelse(mean_mod > 0.35, "white", "black"))
  p_ctx <- ggplot(context_data,
                  aes(x = ctx_up, y = ctx_down, fill = mean_mod)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = n_sites, colour = text_col), size = 2.5) +
    scale_colour_identity() +
    facet_grid(cell_line ~ biotype) +
    scale_fill_gradient(
      low    = "#FFFFFF",
      high   = CONDITION_COLOURS[TREAT_COND],
      name   = "Mean\nmod rate",
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    labs(
      title    = "Trinucleotide context of modified cytosines",
      subtitle = paste0(CONDITION_LABELS[TREAT_COND], " | numbers = site count | NpCpN context"),
      x        = "5' nucleotide (upstream of C)",
      y        = "3' nucleotide (downstream of C)"
    ) +
    theme_srnataps() +
    theme(panel.grid = element_blank())
  save_figure(p_ctx, file.path(opt$figdir, "04e_trinucleotide_context.pdf"),
              width = 10, height = 6)
  message("Figure 4e: Trinucleotide context — done")
}
message("03_modification.R complete.")
