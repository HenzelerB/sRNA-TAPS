# =============================================================================
# 02_biotype.R — Biotype composition figures
#
# Inputs:
#   - 05.biotype_bams/biotype_composition_all_samples.tsv
#
# Outputs:
#   - report/figures/02a_biotype_composition_all.{pdf,png,svg}
#   - report/figures/02b_biotype_composition_mean.{pdf,png,svg}
# =============================================================================

suppressPackageStartupMessages(library(optparse))
source(file.path(Sys.getenv("SRNATAPS_R_DIR", "/mnt/nfs/home/bhenzeler/projects/RNA_TAPS/sRNA-TAPS/srnataps/report/R"), "00_setup.R"))

option_list <- list(
  make_option("--outdir", type = "character"),
  make_option("--figdir", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

biotype_dir  <- file.path(opt$outdir, "05.biotype_bams")
biotype_file <- file.path(biotype_dir, "biotype_composition_all_samples.tsv")

# Assemble from per-sample summary files if combined TSV doesn't exist
if (!file.exists(biotype_file)) {
  summary_files <- list.files(biotype_dir, pattern = "_biotype_summary\\.txt$",
                               full.names = TRUE)
  if (length(summary_files) == 0) stop("No biotype summary files found in ", biotype_dir)
  message("  Assembling biotype composition from ", length(summary_files), " summary files...")
  bio_list <- lapply(summary_files, function(f) {
    sample <- sub("_biotype_summary\\.txt$", "", basename(f))
    df <- read_tsv(f, show_col_types = FALSE) %>%
      dplyr::mutate(
        sample    = sample,
        condition = get_condition(sample),
        cell_line = get_cell_line(sample)
      )
    df
  })
  bio_combined <- dplyr::bind_rows(bio_list)
  write_tsv(bio_combined, biotype_file)
  message("  Written: ", biotype_file)
}

bio <- read_tsv(biotype_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    condition = factor(condition, levels = names(CONDITION_COLOURS)),
    biotype   = factor(biotype,   levels = names(BIOTYPE_COLOURS)),
    cell_line = factor(cell_line, levels = CELL_LINES)
  )


# ── Figure 2a: All samples stacked bar ───────────────────────────────────────

sample_order <- bio %>%
  dplyr::arrange(cell_line, condition, sample) %>%
  dplyr::pull(sample) %>%
  unique()

p_bio_all <- bio %>%
  dplyr::mutate(sample = factor(sample, levels = sample_order)) %>%
  ggplot(aes(x = sample, y = percent, fill = biotype)) +
    geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
    facet_grid(~ cell_line + condition,
               scales = "free_x", space = "free_x",
               labeller = labeller(condition = CONDITION_LABELS)) +
    scale_fill_manual(values = BIOTYPE_COLOURS, name = "RNA biotype") +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(
      title    = "RNA biotype composition",
      subtitle = "Percentage of mapped reads per biotype across all samples",
      x        = NULL,
      y        = "% of mapped reads",
      caption  = "Biotype priority: miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other"
    ) +
    theme_srnataps() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 6.5),
      legend.key.size  = unit(0.4, "cm"),
      strip.text       = element_text(size = 7)
    )

save_figure(p_bio_all, file.path(opt$figdir, "02a_biotype_composition_all.pdf"),
            width = 14, height = 5)
message("Figure 2a: Biotype composition (all samples) — done")


# ── Figure 2b: Mean per condition per cell line ───────────────────────────────

bio_mean <- bio %>%
  dplyr::group_by(biotype, condition, cell_line) %>%
  dplyr::summarise(
    mean_pct = mean(percent, na.rm = TRUE),
    sd_pct   = sd(percent,   na.rm = TRUE),
    n        = dplyr::n(),
    .groups  = "drop"
  )

p_bio_mean <- ggplot(bio_mean,
                     aes(x = condition, y = mean_pct, fill = biotype)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.15,
           position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct),
                position = position_dodge(width = 0.8),
                width = 0.25, colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ cell_line) +
  scale_fill_manual(values = BIOTYPE_COLOURS, name = "RNA biotype") +
  scale_x_discrete(labels = CONDITION_LABELS) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Biotype composition by condition",
    subtitle = "Mean ± SD across biological replicates (n = 3)",
    x        = "Condition",
    y        = "Mean % of mapped reads"
  ) +
  theme_srnataps() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_figure(p_bio_mean, file.path(opt$figdir, "02b_biotype_composition_mean.pdf"),
            width = 8, height = 5)
message("Figure 2b: Biotype composition (mean) — done")
message("02_biotype.R complete.")
