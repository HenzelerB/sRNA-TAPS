# =============================================================================
# 03_modification.R — TAPS modification rate figures
#
# Inputs:
#   - 07.taps_calls/<biotype>/<sample>_<biotype>_taps.tsv
#
# Outputs:
#   - report/figures/03a_modrate_distribution.{pdf,png,svg}
#   - report/figures/03b_top_species.{pdf,png,svg}
#   - report/figures/03c_condition_comparison.{pdf,png,svg}
#   - report/figures/03d_waterfall.{pdf,png,svg}
#   - report/figures/03e_trinucleotide_context.{pdf,png,svg}
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
CALLS_DIR     <- file.path(opt$outdir, "07.taps_calls")

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
  theme_srnataps()

save_figure(p_hist, file.path(opt$figdir, "03a_modrate_distribution.pdf"),
            width = 12, height = 7)

message("Figure 3a: Modification rate distribution — done")

# ── Figure 3b: Top 10 modified sites ─────────────────────────────────────────
# Get top 10 sites per biotype pooled across cell lines
top_site_ids <- calls %>%
  dplyr::filter(condition == TREAT_COND) %>%
  dplyr::mutate(site_id = paste0(chrom, ":", start)) %>%
  dplyr::group_by(biotype, site_id, site_label) %>%
  dplyr::summarise(median_mod = median(mod_rate), .groups = "drop") %>%
  dplyr::group_by(biotype) %>%
  dplyr::slice_max(median_mod, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  # Make labels unique at this stage
  dplyr::group_by(biotype, site_label) %>%
  dplyr::mutate(
    n_pos = dplyr::n_distinct(site_id),
    site_label_uniq = ifelse(n_pos > 1,
                             paste0(site_label, " (", site_id, ")"),
                             site_label)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(biotype, site_id, site_label, site_label_uniq)

# Get per-cell-line stats for those top 10 sites
top_sites <- calls %>%
  dplyr::filter(condition == TREAT_COND) %>%
  dplyr::mutate(site_id = paste0(chrom, ":", start)) %>%
  dplyr::inner_join(top_site_ids, by = c("biotype", "site_id", "site_label")) %>%
  dplyr::group_by(biotype, site_id, site_label, site_label_uniq, cell_line) %>%
  dplyr::summarise(
    median_mod = median(mod_rate),
    mean_cov   = mean(coverage),
    .groups    = "drop"
  )

if (nrow(top_sites) > 0) {
  # 2x2 matrix layout: miRNA+tRNA on left, rRNA+snoRNA on right
  top_matrix <- top_sites %>%
    dplyr::filter(biotype %in% c("miRNA", "tRNA", "rRNA", "snoRNA")) %>%
    dplyr::mutate(
      col_group = ifelse(biotype %in% c("miRNA", "tRNA"), "miRNA & tRNA", "rRNA & snoRNA"),
      col_group = factor(col_group, levels = c("miRNA & tRNA", "rRNA & snoRNA")),
      biotype   = factor(biotype, levels = c("miRNA", "tRNA", "rRNA", "snoRNA")),
      site_label_uniq = reorder(site_label_uniq, median_mod)
    )

  p_top <- ggplot(top_matrix,
                  aes(x = median_mod * 100, y = site_label_uniq,
                      size = log10(mean_cov + 1),
                      colour = biotype, shape = cell_line)) +
    geom_segment(aes(x = 0, xend = median_mod * 100,
                     y = site_label_uniq, yend = site_label_uniq),
                 colour = "grey85", linewidth = 0.3) +
    geom_point(alpha = 0.85) +
    facet_wrap(~ biotype, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_shape_manual(values = c("Caco2" = 16, "HEK" = 17), name = "Cell line") +
    scale_size_continuous(name = "log10(coverage)", range = c(1, 5)) +
    scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
    labs(
      title    = paste0("Top 10 modified sites per biotype (", CONDITION_LABELS[TREAT_COND], ")"),
      subtitle = "Caco2 (circle) and HEK (triangle) — ranked by median modification rate",
      x        = "Modification rate (%)",
      y        = NULL
    ) +
    theme_srnataps() +
    theme(
      axis.text.y     = element_text(size = 7),
      legend.position = "bottom"
    )

  save_figure(p_top, file.path(opt$figdir, "03b_top_sites.pdf"),
              width = 8, height = 4.8)
  message("Figure 3b: Top modified sites — done")
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
  cond_wide <- cond_wide %>%
    dplyr::mutate(above_diag = .data[[TREAT_COND]] > .data[[CTRL_COND]])

  p_scatter <- ggplot(cond_wide,
                      aes(x = .data[[CTRL_COND]], y = .data[[TREAT_COND]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.4) +
    geom_point(data = dplyr::filter(cond_wide, !above_diag),
               colour = "grey80", alpha = 0.25, size = 0.8, shape = 16) +
    geom_point(data = dplyr::filter(cond_wide,  above_diag),
               aes(colour = cell_line),
               alpha = 0.35, size = 1.1, shape = 16) +
    facet_wrap(~ biotype, nrow = 2) +
    scale_colour_manual(values = CELL_COLOURS, name = "Cell line") +
    scale_x_continuous(labels = function(x) x * 100, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(
      title    = "Condition comparison: PB-only vs TET+PB",
      subtitle = "Caco2 (mint) and HEK (gold) | above diagonal = TET-enriched",
      x        = paste0(CONDITION_LABELS[CTRL_COND],  " mod rate (%)"),
      y        = paste0(CONDITION_LABELS[TREAT_COND], " mod rate (%)"),
      caption  = "Dashed line = identity"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    theme_srnataps() +
    theme(legend.position = "bottom")
  save_figure(p_scatter, file.path(opt$figdir, "03c_condition_comparison.pdf"),
              width = 7, height = 5)
  message("Figure 3c: Condition comparison scatter — done")
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
  save_figure(p_waterfall, file.path(opt$figdir, "03d_waterfall.pdf"),
              width = 9, height = 5)
  message("Figure 3d: Waterfall plot — done")
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
  save_figure(p_ctx, file.path(opt$figdir, "03e_trinucleotide_context.pdf"),
              width = 10, height = 6)
  message("Figure 3e: Trinucleotide context — done")
}
message("03_modification.R complete.")
