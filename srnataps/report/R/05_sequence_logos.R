# =============================================================================
# 05_sequence_logos.R — Sequence logos around modified cytosines
#
# Extracts genomic sequence windows around high-confidence TAPS m5C sites
# and generates sequence logos showing nucleotide context preference.
#
# Requires Bioconductor packages:
#   BiocManager::install(c("Rsamtools", "GenomicRanges"))
#   install.packages("ggseqlogo")
#
# Inputs:
#   - 07.taps_calls/<biotype>/<sample>_<biotype>_taps.tsv
#     Required columns: chrom, start, strand, mod_rate, padj, snp_flag
#
# Outputs (per window size ±5 and ±10, per biotype + combined):
#   - report/figures/05_seqlogo_combined_w5.{pdf,png,svg}
#   - report/figures/05_seqlogo_combined_w10.{pdf,png,svg}
#   - report/figures/05_seqlogo_<biotype>_w5.{pdf,png,svg}
#   - report/figures/05_seqlogo_<biotype>_w10.{pdf,png,svg}
#
# Site selection (high confidence):
#   padj    < 0.01   (BH-corrected binomial test)
#   mod_rate > 0.50   (at least 50% C→T conversion)
#   snp_flag == PASS  (not a SNP)
#   coverage >= min_cov
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggseqlogo)
  library(Rsamtools)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(stringr)
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
  make_option("--outdir",    type = "character", help = "Project output directory"),
  make_option("--figdir",    type = "character", default = NULL),
  make_option("--calls-dir", type = "character", default = NULL),
  make_option("--genome", type = "character",
              help = "Indexed reference FASTA used by the pipeline"),
  make_option("--biotypes",  type = "character",
              default = "miRNA,tRNA,rRNA,snoRNA",
              help = "Comma-separated biotypes [default: miRNA,tRNA,rRNA,snoRNA]"),
  make_option("--min-cov",   type = "integer",   default = 5),
  make_option("--min-mod",   type = "double",    default = 0.50,
              help = "Minimum mod_rate for high-confidence sites [default: 0.50]"),
  make_option("--max-padj",  type = "double",    default = 0.01,
              help = "Maximum adjusted p-value [default: 0.01]"),
  make_option("--max-sites", type = "integer",   default = 2000,
              help = "Max sites per biotype (random sample if more) [default: 2000]"),
  make_option("--condition", type = "character", default = "treat",
              help = "Condition to use [default: treat]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

BIOTYPES_USE <- strsplit(opt$biotypes, ",")[[1]]
CALLS_DIR    <- if (is.null(opt$`calls-dir`)) {
  file.path(opt$outdir, "07.taps_calls")
} else opt$`calls-dir`
if (is.null(opt$genome) || !file.exists(opt$genome)) {
  stop("--genome must point to the indexed pipeline reference FASTA")
}
GENOME       <- FaFile(opt$genome)
AVAILABLE_SEQS <- as.character(seqnames(scanFaIndex(GENOME)))
WINDOWS      <- c(5, 10)   # ±5 nt and ±10 nt


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load high-confidence sites
# ══════════════════════════════════════════════════════════════════════════════

load_hc_sites <- function(calls_dir, biotypes, condition, min_cov, min_mod, max_padj) {
  message("Loading high-confidence sites...")
  message("  Filters: mod_rate >= ", min_mod,
          " | padj < ", max_padj,
          " | coverage >= ", min_cov,
          " | snp_flag == PASS")

  all_sites <- lapply(biotypes, function(bt) {
    bt_dir <- file.path(calls_dir, bt)
    if (!dir.exists(bt_dir)) return(NULL)

    files <- list.files(bt_dir, pattern = "_taps\\.tsv$", full.names = TRUE)
    # Filter to target condition, exclude untreated and PB-only aliases.
    files <- files[grepl(condition, files) &
                   !grepl("no[-_]?treat|untreated|untr", files, ignore.case = TRUE) &
                   !grepl("pb[-_]?ctrl|pb[-_]?control|pb[-_]?only", files, ignore.case = TRUE)]

    if (length(files) == 0) return(NULL)

    lapply(files, function(f) {
      sample <- sub(paste0("_", bt, "_taps\\.tsv$"), "", basename(f))
      tryCatch({
        df <- read_tsv(f, show_col_types = FALSE,
                             col_types = cols(chrom = col_character())) %>%
          dplyr::filter(
            snp_flag == "PASS",
            coverage >= min_cov,
            mod_rate >= min_mod,
            padj     <  max_padj
          ) %>%
          dplyr::mutate(
            biotype   = bt,
            sample    = sample,
            cell_line = get_cell_line(sample)
          )
        df
      }, error = function(e) NULL)
    }) %>% dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows(Filter(Negate(is.null), .))

  message("  Total high-confidence sites: ", nrow(all_sites))
  all_sites
}

sites <- load_hc_sites(CALLS_DIR, BIOTYPES_USE, opt$condition,
                       opt$`min-cov`, opt$`min-mod`, opt$`max-padj`)

if (nrow(sites) == 0) {
  stop("No high-confidence sites found. Check filters or that calling has completed.")
}

# Print site counts per biotype
site_counts <- sites %>%
  dplyr::count(biotype, cell_line) %>%
  dplyr::arrange(biotype, cell_line)
message("\nSites per biotype:")
print(as.data.frame(site_counts))


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Extract genomic sequences around each site
# ══════════════════════════════════════════════════════════════════════════════

# Match Ensembl/UCSC chromosome naming to the configured FASTA.
normalise_chrom <- function(chrom) {
  vapply(chrom, function(value) {
    candidates <- c(
      value,
      if (grepl("^chr", value)) sub("^chr", "", value) else paste0("chr", value),
      if (value %in% c("M", "MT", "chrM", "chrMT")) {
        c("MT", "M", "chrM", "chrMT")
      }
    )
    match <- candidates[candidates %in% AVAILABLE_SEQS]
    if (length(match)) match[[1]] else value
  }, character(1))
}

extract_sequences <- function(sites_df, genome, half_window) {
  # half_window = 5 means we extract from pos-5 to pos+5 (11 nt total)
  # The modified C is at position half_window + 1 (1-based, centre)

  sites_df <- sites_df %>%
    dplyr::mutate(
      chrom_norm = normalise_chrom(chrom),
      seq_start  = start - half_window,     # 0-based start
      seq_end    = start + half_window + 1  # 0-based end (exclusive)
    ) %>%
    dplyr::filter(
      seq_start >= 0,
      chrom_norm %in% seqnames(genome)
    )

  if (nrow(sites_df) == 0) return(character(0))

  # Limit to max_sites per call to avoid memory issues
  if (nrow(sites_df) > opt$`max-sites`) {
    set.seed(42)
    sites_df <- dplyr::slice_sample(sites_df, n = opt$`max-sites`)
    message("  Subsampled to ", opt$`max-sites`, " sites")
  }

  message("  Extracting sequences for ", nrow(sites_df),
          " sites (window ±", half_window, " nt)...")

  # Build GRanges and fetch sequences
  # Handle missing strand column gracefully
  strand_vec <- if ("strand" %in% names(sites_df)) {
    ifelse(!is.na(sites_df$strand) & sites_df$strand == "-", "-", "+")
  } else {
    rep("+", nrow(sites_df))
  }

  gr <- GRanges(
    seqnames = sites_df$chrom_norm,
    ranges   = IRanges(start = sites_df$seq_start + 1,
                       end   = sites_df$seq_start + (2 * half_window + 1)),
    strand   = strand_vec
  )

  seqs <- tryCatch(
    as.character(getSeq(genome, gr)),
    error = function(e) {
      message("  WARNING: getSeq failed: ", e$message)
      character(0)
    }
  )

  # Keep only sequences of correct length (some near chromosome ends may be shorter)
  expected_len <- 2 * half_window + 1
  seqs <- seqs[nchar(seqs) == expected_len]

  # Remove sequences with N
  seqs <- seqs[!grepl("N", seqs, ignore.case = TRUE)]
  seqs <- toupper(seqs)

  message("  Valid sequences retained: ", length(seqs))
  seqs
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Build sequence logo per window per biotype
# ══════════════════════════════════════════════════════════════════════════════

# Custom DNA colour scheme matching TAPS context
# C is highlighted at the centre — the modified base
dna_colours <- make_col_scheme(
  chars = c("A", "C", "G", "T"),
  cols  = c("#2ECC71", "#E74C3C", "#F39C12", "#3498DB")
)

make_logo_plot <- function(seqs, title, half_window, n_sites) {
  if (length(seqs) < 10) {
    message("  Too few sequences (", length(seqs), ") for logo — skipping")
    return(NULL)
  }

  total_width <- 2 * half_window + 1
  centre_pos  <- half_window + 1

  # x-axis labels: -hw ... -1, C, +1 ... +hw
  x_labels <- c(
    paste0("-", half_window:1),
    "C",
    paste0("+", 1:half_window)
  )

  p <- ggplot() +
    geom_logo(seqs,
              method    = "bits",
              col_scheme = dna_colours,
              font      = "roboto_medium") +
    annotate("rect",
             xmin = centre_pos - 0.5,
             xmax = centre_pos + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "#E74C3C", alpha = 0.08) +
    annotate("segment",
             x = centre_pos, xend = centre_pos,
             y = 0, yend = 2,
             colour = "#E74C3C", linewidth = 0.4, linetype = "dashed") +
    scale_x_continuous(
      breaks = 1:total_width,
      labels = x_labels,
      expand = c(0.01, 0.01)
    ) +
    labs(
      title    = title,
      subtitle = paste0("n = ", scales::comma(n_sites),
                        " high-confidence sites (mod_rate ≥ ",
                        opt$`min-mod`, ", padj < ", opt$`max-padj`, ")"),
      x        = paste0("Position relative to m5C (±", half_window, " nt)"),
      y        = "Bits"
    ) +
    theme_arial8 +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 6, family = "Arial")
    )

  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Generate and save all logos
# ══════════════════════════════════════════════════════════════════════════════

for (hw in WINDOWS) {
  message("\n", strrep("=", 55))
  message("Window ±", hw, " nt")
  message(strrep("=", 55))

  # ── Combined logo (all biotypes) ───────────────────────────────────────────
  message("Combined (all biotypes):")
  seqs_all <- extract_sequences(sites, GENOME, hw)

  p_combined <- make_logo_plot(
    seqs  = seqs_all,
    title = paste0("Sequence context of m5C sites — all biotypes (±", hw, " nt)"),
    half_window = hw,
    n_sites     = length(seqs_all)
  )

  if (!is.null(p_combined)) {
    save_figure(
      p_combined,
      file.path(opt$figdir, paste0("05_seqlogo_combined_w", hw, ".pdf")),
      width = ifelse(hw == 5, 5, 8),
      height = 2.5
    )
    message("  Saved: 05_seqlogo_combined_w", hw)
  }

  # ── Per-biotype logos ─────────────────────────────────────────────────────
  bio_plots <- list()

  for (bt in BIOTYPES_USE) {
    bt_sites <- sites %>% dplyr::filter(biotype == bt)
    if (nrow(bt_sites) < 10) {
      message("  ", bt, ": too few sites (", nrow(bt_sites), ") — skipping")
      next
    }

    message(bt, " (", nrow(bt_sites), " sites):")
    seqs_bt <- extract_sequences(bt_sites, GENOME, hw)

    p_bt <- make_logo_plot(
      seqs        = seqs_bt,
      title       = paste0(bt, " m5C context (±", hw, " nt)"),
      half_window = hw,
      n_sites     = length(seqs_bt)
    )

    if (!is.null(p_bt)) {
      # Individual figure
      save_figure(
        p_bt,
        file.path(opt$figdir,
                  paste0("05_seqlogo_", bt, "_w", hw, ".pdf")),
        width  = ifelse(hw == 5, 5, 8),
        height = 2.5
      )
      bio_plots[[bt]] <- p_bt
    }
  }

  # ── Panel figure: all biotypes in one PDF ─────────────────────────────────
  if (length(bio_plots) >= 2) {
    panel <- wrap_plots(bio_plots, ncol = 1) +
      plot_annotation(
        title    = paste0("m5C sequence context per RNA biotype (±", hw, " nt)"),
        subtitle = paste0("High-confidence sites: mod_rate ≥ ", opt$`min-mod`,
                          " | padj < ", opt$`max-padj`,
                          " | TET+PB condition"),
        theme    = theme_arial8
      )

    save_figure(
      panel,
      file.path(opt$figdir,
                paste0("05_seqlogo_per_biotype_panel_w", hw, ".pdf")),
      width  = ifelse(hw == 5, 5, 8),
      height = 2.5 * length(bio_plots)
    )
    message("Panel figure saved: 05_seqlogo_per_biotype_panel_w", hw)
  }
}

# ── Summary table of sites used ───────────────────────────────────────────────
summary_out <- file.path(opt$figdir, "05_seqlogo_site_summary.tsv")
sites %>%
  dplyr::count(biotype, cell_line, name = "n_sites") %>%
  dplyr::arrange(biotype, cell_line) %>%
  readr::write_tsv(summary_out)
message("\nSite summary: ", summary_out)
message("05_sequence_logos.R complete.")
