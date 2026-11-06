#!/usr/bin/env Rscript
# =============================================================================
# 07_global_distribution.R — Figure 3: global modification-rate percentiles
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
  make_option("--calls-dir", type = "character", default = NULL),
  make_option("--min-cov", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$outdir)) stop("--outdir is required")
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(opt$`calls-dir`)) {
  opt$`calls-dir` <- file.path(opt$outdir, "07.taps_calls_3letter")
}
call_files <- list.files(
  opt$`calls-dir`, pattern = "_taps\\.tsv$",
  recursive = TRUE, full.names = TRUE
)
if (length(call_files) == 0) stop("No stage-07 call files found")

condition_order <- unique(c(UNTR_COND, CTRL_COND, TREAT_COND))
biotype_order <- c("rRNA", "miRNA", "snoRNA", "tRNA", "snRNA",
                   "piRNA", "lncRNA", "other")

parse_file <- function(f) {
  biotype <- basename(dirname(f))
  bn <- basename(f)
  sample <- sub(paste0("_", biotype, "_taps\\.tsv$"), "", bn)

  tryCatch({
    d <- readr::read_tsv(f, show_col_types = FALSE)
    required_columns <- c("coverage", "mod_rate", "snp_flag")
    if (!all(required_columns %in% names(d))) return(NULL)
    d %>%
      dplyr::mutate(
        condition = get_condition(sample),
        cell_line = get_cell_line(sample),
        biotype = biotype
      ) %>%
      dplyr::filter(coverage >= opt$`min-cov`, snp_flag == "PASS", is.finite(mod_rate),
                    mod_rate >= 0, mod_rate <= 1)
  }, error = function(e) NULL)
}

rates <- lapply(call_files, parse_file) %>%
  Filter(Negate(is.null), .) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(condition %in% condition_order,
                biotype %in% biotype_order) %>%
  dplyr::mutate(
    condition = factor(condition, levels = condition_order),
    biotype = factor(biotype, levels = biotype_order)
  )
if (nrow(rates) == 0) stop("No eligible stringent sites after filtering")

percentiles_observed <- rates %>%
  dplyr::group_by(biotype, condition) %>%
  dplyr::summarise(
    n_sites = dplyr::n(),
    q1 = stats::quantile(mod_rate, 0.25, na.rm = TRUE),
    median = stats::median(mod_rate, na.rm = TRUE),
    q3 = stats::quantile(mod_rate, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

percentiles <- tidyr::expand_grid(
  biotype = factor(biotype_order, levels = biotype_order),
  condition = factor(condition_order, levels = condition_order)
) %>%
  dplyr::left_join(percentiles_observed, by = c("biotype", "condition")) %>%
  dplyr::mutate(
    annotation = dplyr::if_else(
      is.na(n_sites), "No eligible\nstage-07 sites",
      sprintf("Md: %.1f%%\nQ1: %.1f%%\nQ3: %.1f%%",
              median * 100, q1 * 100, q3 * 100)
    ),
    annotation_y = dplyr::if_else(is.na(n_sites), 0.52, 1.07)
  )

readr::write_tsv(
  percentiles,
  file.path(opt$figdir, "03_global_modrate_percentiles.tsv")
)

p_global <- ggplot(rates, aes(x = condition, y = mod_rate,
                              fill = condition, colour = condition)) +
  geom_violin(
    trim = TRUE, scale = "width", width = 0.88,
    alpha = 0.72, linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.12, outlier.shape = NA,
    fill = "white", colour = "#1C4062", linewidth = 0.35
  ) +
  geom_text(
    data = percentiles,
    aes(x = condition, y = annotation_y, label = annotation),
    inherit.aes = FALSE,
    size = 1.85, lineheight = 0.88, colour = "#1C4062",
    fontface = "plain"
  ) +
  facet_wrap(~ biotype, ncol = 4, drop = FALSE) +
  scale_fill_manual(values = CONDITION_COLOURS,
                    labels = CONDITION_LABELS, name = "Condition") +
  scale_colour_manual(values = CONDITION_COLOURS,
                      labels = CONDITION_LABELS, guide = "none") +
  scale_x_discrete(labels = CONDITION_LABELS) +
  scale_y_continuous(
    limits = c(0, 1.17),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Global modification-rate distribution by RNA biotype",
    subtitle = paste0(
      "Stage-07 PASS sites with coverage ≥ ", opt$`min-cov`,
      "× | Violin density with median and interquartile range"
    ),
    x = NULL, y = "Modification rate",
    caption = "Md = median | Q1 = 25th percentile | Q3 = 75th percentile"
  ) +
  theme_srnataps(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "#1C4062", colour = NA),
    strip.text = element_text(colour = "white", face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 6.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

save_figure(
  p_global,
  file.path(opt$figdir, "03_global_modrate_distribution.pdf"),
  width = 12, height = 7.5
)

message("Figure 3: Global modification-rate distributions — ",
  nrow(rates), " stage-07 condition measurements across ",
        dplyr::n_distinct(rates$biotype), " RNA biotypes — done")
