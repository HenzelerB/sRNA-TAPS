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
suppressPackageStartupMessages(library(ggdist))
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


# ── Figure 2: Biotype composition — one panel per biotype ───────────────────

bio_mean <- bio %>%
  dplyr::group_by(biotype, condition, cell_line) %>%
  dplyr::summarise(
    mean_pct = mean(percent, na.rm = TRUE),
    sd_pct   = sd(percent,   na.rm = TRUE),
    .groups  = "drop"
  )

bio_indiv <- bio %>%
  dplyr::mutate(condition = factor(condition, levels = names(CONDITION_COLOURS)))

p_bio <- ggplot(bio_mean,
               aes(x = condition, y = mean_pct, fill = biotype)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  facet_wrap(~ cell_line) +
  scale_fill_manual(values = BIOTYPE_COLOURS, name = "RNA biotype") +
  scale_x_discrete(labels = CONDITION_LABELS) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "RNA biotype composition by condition",
    subtitle = "Mean percentage of mapped reads per biotype",
    x        = NULL,
    y        = "% of mapped reads",
    caption  = "Biotype priority: miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other"
  ) +
  theme_srnataps() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 30, hjust = 1)
  )

save_figure(p_bio, file.path(opt$figdir, "02_biotype_composition.pdf"),
            width = 5, height = 5)
message("Figure 2: Biotype composition — done")
message("02_biotype.R complete.")
