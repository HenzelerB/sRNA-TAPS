# =============================================================================
# 04_benchmark.R — Benchmarking comparison figures
#
# Inputs:
#   - 09.compare/concordance_summary.tsv
#   - 09.compare/correlation_summary.tsv
#   - 09.compare/shared_*_*.tsv
#
# Outputs:
#   - report/figures/04a_concordance_heatmap.{pdf,png,svg}
#   - report/figures/04b_correlation_scatter.{pdf,png,svg}
#   - report/figures/04c_venn.{pdf,png,svg}
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})
source(file.path(Sys.getenv("SRNATAPS_R_DIR", "/mnt/nfs/home/bhenzeler/projects/RNA_TAPS/sRNA-TAPS/srnataps/report/R"), "00_setup.R"))

option_list <- list(
  make_option("--outdir", type = "character"),
  make_option("--figdir", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

COMPARE_DIR <- file.path(opt$outdir, "09.compare")

TOOL_LABELS <- c(
  "rastair"  = "rastair\n(Bowtie1, CpG)",
  "astair"   = "asTair\n(Bowtie1, all ctx)",
  "bismark"  = "Bismark\n(bowtie2, inverted)"
)

TOOL_COLOURS <- c(
  "rastair"  = "#2196F3",
  "astair"   = "#4CAF50",
  "bismark"  = "#FF9800"
)


# ── Figure 4a: Concordance heatmap ───────────────────────────────────────────

conc_file <- file.path(COMPARE_DIR, "concordance_summary.tsv")
if (file.exists(conc_file)) {
  conc <- read_tsv(conc_file, show_col_types = FALSE) %>%
    dplyr::filter(grepl("treat", condition) & !grepl("no.treat|no_treat", condition)) %>%
    dplyr::mutate(
      tool    = factor(tool,    levels = names(TOOL_LABELS)),
      biotype = factor(biotype, levels = names(BIOTYPE_COLOURS))
    )

  p_conc <- ggplot(conc, aes(x = tool, y = biotype, fill = jaccard)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", jaccard)),
              size = 3, colour = "grey10") +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         name = "Jaccard\nindex",
                         limits = c(0, 1)) +
    scale_x_discrete(labels = TOOL_LABELS) +
    labs(
      title    = "Site-level concordance: custom pipeline vs benchmarks",
      subtitle = "TET+PB condition | Jaccard index at called sites",
      x        = NULL,
      y        = "RNA biotype",
      caption  = "Bismark chemistry is inverted (1 - rate) before comparison"
    ) +
    theme_srnataps() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 8))

  save_figure(p_conc, file.path(opt$figdir, "04a_concordance_heatmap.pdf"),
              width = 7, height = 6)
  message("Figure 4a: Concordance heatmap — done")
} else {
  message("WARNING: concordance_summary.tsv not found — skipping Figure 4a")
}


# ── Figure 4b: Correlation scatter (mod_rate at shared sites) ────────────────

corr_file <- file.path(COMPARE_DIR, "correlation_summary.tsv")
if (file.exists(corr_file)) {
  corr <- read_tsv(corr_file, show_col_types = FALSE) %>%
    dplyr::filter(condition == "treat", !is.na(pearson_r)) %>%
    dplyr::mutate(tool = factor(tool, levels = names(TOOL_LABELS)))

  p_corr <- ggplot(corr, aes(x = biotype, y = pearson_r,
                              colour = tool, group = tool)) +
    geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_point(aes(shape = tool), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(values = TOOL_COLOURS, labels = TOOL_LABELS,
                        name = "Tool") +
    scale_shape_manual(values = c(16, 17, 15), labels = TOOL_LABELS,
                       name = "Tool") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(
      title    = "Pearson correlation of modification rates at shared sites",
      subtitle = "TET+PB condition per biotype",
      x        = "RNA biotype",
      y        = "Pearson r",
      caption  = "Negative correlation for Bismark is expected (chemistry inversion)"
    ) +
    theme_srnataps() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_figure(p_corr, file.path(opt$figdir, "04b_correlation.pdf"),
              width = 9, height = 5)
  message("Figure 4b: Correlation plot — done")

  # Also make scatter plots for each tool at the treat condition, miRNA biotype
  shared_files <- list.files(COMPARE_DIR,
                              pattern = "shared_treat.*\\.tsv$",
                              full.names = TRUE)
  if (length(shared_files) > 0) {
    shared_data <- lapply(shared_files, function(f) {
      tool    <- sub(".*shared_treat_[^_]+_([^.]+)\\.tsv", "\\1", basename(f))
      biotype <- sub(".*shared_treat_([^_]+)_.*\\.tsv", "\\1", basename(f))
      df <- read_tsv(f, show_col_types = FALSE,
                     col_types = cols(site_key = col_character()))
      df$tool    <- tool
      df$biotype <- biotype
      df
    }) %>% dplyr::bind_rows()

    if (nrow(shared_data) > 0 &&
        "mod_rate_custom" %in% names(shared_data)) {
      tool_col <- grep("mod_rate_(?!custom)", names(shared_data),
                       perl = TRUE, value = TRUE)[1]
      if (!is.na(tool_col)) {
        p_shared <- ggplot(shared_data,
                           aes_string(x = "mod_rate_custom",
                                      y = tool_col,
                                      colour = "biotype")) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                      colour = "grey60", linewidth = 0.3) +
          geom_point(alpha = 0.4, size = 0.8) +
          facet_wrap(~ tool + biotype, ncol = 4) +
          scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
          scale_x_continuous(labels = percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          scale_y_continuous(labels = percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          labs(
            title = "Modification rate at shared sites",
            x     = "Custom pipeline mod rate",
            y     = "Benchmark tool mod rate"
          ) +
          theme_srnataps(base_size = 9)

        save_figure(p_shared,
                    file.path(opt$figdir, "04b_shared_sites_scatter.pdf"),
                    width = 14, height = 8)
        message("Figure 4b extra: Shared sites scatter — done")
      }
    }
  }
}


# ── Figure 4c: n_shared bar chart (Venn-style summary) ───────────────────────

if (file.exists(conc_file)) {
  conc_bar <- read_tsv(conc_file, show_col_types = FALSE) %>%
    dplyr::filter(grepl("treat", condition) & !grepl("no.treat|no_treat", condition)) %>%
    dplyr::select(tool, biotype, sites_custom, sites_tool, shared) %>%
    tidyr::pivot_longer(cols = c(sites_custom, sites_tool, shared),
                        names_to = "category", values_to = "n_sites") %>%
    dplyr::mutate(
      category = dplyr::recode(category,
        "sites_custom" = "Custom only",
        "sites_tool"   = "Tool only",
        "shared"       = "Shared"
      ),
      category = factor(category,
                        levels = c("Custom only", "Shared", "Tool only")),
      tool     = factor(tool, levels = names(TOOL_LABELS))
    )

  if (nrow(conc_bar) == 0 || length(unique(conc_bar$biotype)) == 0) {
    message("WARNING: No data for Figure 4c — skipping")
  } else {
  p_venn <- ggplot(conc_bar, aes(x = tool, y = n_sites, fill = category)) +
    geom_col(position = "dodge", width = 0.75,
             colour = "grey30", linewidth = 0.2) +
    facet_wrap(~ biotype, scales = "free_y") +
    scale_fill_manual(
      values = c("Custom only" = "#E41A1C",
                 "Shared"      = "#4DAF4A",
                 "Tool only"   = "#377EB8"),
      name = "Site category"
    ) +
    scale_x_discrete(labels = c("rastair" = "rastair",
                                 "astair"  = "asTair",
                                 "bismark" = "Bismark")) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      title    = "Called sites: custom pipeline vs benchmark tools",
      subtitle = "TET+PB condition per biotype",
      x        = NULL,
      y        = "Number of sites"
    ) +
    theme_srnataps() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

  save_figure(p_venn, file.path(opt$figdir, "04c_site_overlap.pdf"),
              width = 12, height = 8)
  message("Figure 4c: Site overlap bar chart — done")
  } # end guard
}

message("04_benchmark.R complete.")
