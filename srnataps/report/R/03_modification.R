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
source(file.path(Sys.getenv("SRNATAPS_R_DIR", "/mnt/nfs/home/bhenzeler/projects/RNA_TAPS/sRNA-TAPS/srnataps/report/R"), "00_setup.R"))

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
    list.files(file.path(calls_dir, bt), pattern = "_taps\\.tsv$",
               full.names = TRUE)
  }))

  if (length(all_files) == 0) stop("No TAPS call files found in ", calls_dir)
  message("  Loading ", length(all_files), " call files...")

  lapply(all_files, function(f) {
    bt     <- basename(dirname(f))
    sample <- sub(paste0("_", bt, "_taps\\.tsv$"), "", basename(f))
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE) %>%
        dplyr::filter(coverage >= min_cov, snp_flag == "PASS") %>%
        dplyr::mutate(
          biotype   = bt,
          sample    = sample,
          condition = get_condition(sample),
          cell_line = get_cell_line(sample)
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


# ── Figure 3a: Modification rate distribution (violin + boxplot) ──────────────

p_dist <- calls %>%
  dplyr::filter(condition %in% c("no_treat", "pb_ctrl", "treat")) %>%
  ggplot(aes(x = condition, y = mod_rate, fill = condition)) +
    geom_violin(trim = TRUE, alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.4, colour = "grey20") +
    facet_wrap(~ biotype, nrow = 1) +
    scale_fill_manual(values = CONDITION_COLOURS, labels = CONDITION_LABELS,
                      name = "Condition") +
    scale_x_discrete(labels = CONDITION_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(
      title    = "TAPS modification rate distribution",
      subtitle = "Per-site mod_rate at PASS positions (SNP-filtered, BH-corrected)",
      x        = NULL,
      y        = "Modification rate",
      caption  = paste0("Min coverage: ", opt$`min-cov`, "x | SNP_flag == PASS only")
    ) +
    theme_srnataps() +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )

save_figure(p_dist, file.path(opt$figdir, "03a_modrate_distribution.pdf"),
            width = 12, height = 5)
message("Figure 3a: Modification rate distribution — done")


# ── Figure 3b: Top 20 modified RNA species (dot plot) ────────────────────────

top_species <- calls %>%
  dplyr::filter(condition == "treat", !is.na(gene_id), gene_id != ".") %>%
  dplyr::group_by(biotype, gene_id, cell_line) %>%
  dplyr::summarise(
    median_mod = median(mod_rate),
    mean_cov   = mean(coverage),
    n_sites    = n(),
    .groups    = "drop"
  ) %>%
  dplyr::filter(n_sites >= 2) %>%
  dplyr::group_by(biotype, cell_line) %>%
  dplyr::slice_max(median_mod, n = 20) %>%
  dplyr::ungroup()

if (nrow(top_species) > 0) {
  p_top <- top_species %>%
    dplyr::mutate(gene_id = reorder(gene_id, median_mod)) %>%
    ggplot(aes(x = median_mod, y = gene_id,
               size = log10(mean_cov + 1), colour = biotype)) +
      geom_point(alpha = 0.8) +
      facet_grid(biotype ~ cell_line, scales = "free_y", space = "free_y") +
      scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
      scale_size_continuous(name = "log10(mean coverage)",
                            range = c(1, 5), breaks = c(1, 2, 3)) +
      scale_x_continuous(labels = percent_format(accuracy = 1),
                         limits = c(0, 1)) +
      labs(
        title    = "Top modified RNA species (TET+PB condition)",
        subtitle = "Top 20 per biotype ranked by median modification rate",
        x        = "Median modification rate",
        y        = NULL
      ) +
      theme_srnataps() +
      theme(axis.text.y = element_text(size = 6))

  save_figure(p_top, file.path(opt$figdir, "03b_top_species.pdf"),
              width = 10, height = 12)
  message("Figure 3b: Top modified species — done")
}


# ── Figure 3c: Condition comparison scatter (PB vs TET+PB) ───────────────────

cond_wide <- calls %>%
  dplyr::filter(condition %in% c("pb_ctrl", "treat"),
                !is.na(gene_id), gene_id != ".") %>%
  dplyr::group_by(biotype, gene_id, cell_line, condition) %>%
  dplyr::summarise(median_mod = median(mod_rate), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = median_mod) %>%
  dplyr::filter(!is.na(pb_ctrl), !is.na(treat))

if (nrow(cond_wide) > 0) {
  p_scatter <- ggplot(cond_wide, aes(x = pb_ctrl, y = treat, colour = biotype)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(
      data = dplyr::filter(cond_wide, treat > 0.3),
      aes(label = gene_id),
      size = 2.5, max.overlaps = 15, segment.colour = "grey70"
    ) +
    facet_grid(biotype ~ cell_line) +
    scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title    = "Condition comparison: PB-only vs TET+PB",
      subtitle = "Sites above diagonal = enriched by TET oxidation (genuine m5C/5hmC)",
      x        = "PB-only median mod rate",
      y        = "TET+PB median mod rate",
      caption  = "Dashed line: y = x (no enrichment)"
    ) +
    theme_srnataps()

  save_figure(p_scatter, file.path(opt$figdir, "03c_condition_comparison.pdf"),
              width = 8, height = 10)
  message("Figure 3c: Condition comparison scatter — done")
}


# ── Figure 3d: Waterfall plot (sites ranked by mod_rate) ─────────────────────

waterfall_data <- calls %>%
  dplyr::filter(condition == "treat") %>%
  dplyr::arrange(desc(mod_rate)) %>%
  dplyr::mutate(rank = row_number())

if (nrow(waterfall_data) > 0) {
  p_waterfall <- waterfall_data %>%
    dplyr::slice_head(n = 500) %>%   # top 500 sites for clarity
    ggplot(aes(x = rank, y = mod_rate, colour = biotype)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_colour_manual(values = BIOTYPE_COLOURS, name = "Biotype") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title    = "Top 500 modification sites (TET+PB)",
        subtitle = "Ranked by modification rate, coloured by RNA biotype",
        x        = "Rank",
        y        = "Modification rate"
      ) +
      theme_srnataps()

  save_figure(p_waterfall, file.path(opt$figdir, "03d_waterfall.pdf"),
              width = 9, height = 5)
  message("Figure 3d: Waterfall plot — done")
}


# ── Figure 3e: Trinucleotide context heatmap ─────────────────────────────────

context_data <- calls %>%
  dplyr::filter(condition == "treat", nchar(context) == 3) %>%
  dplyr::mutate(
    ctx_up   = substr(context, 1, 1),
    ctx_down = substr(context, 3, 3)
  ) %>%
  dplyr::group_by(biotype, ctx_up, ctx_down) %>%
  dplyr::summarise(
    mean_mod  = mean(mod_rate),
    n_sites   = n(),
    .groups   = "drop"
  ) %>%
  dplyr::filter(n_sites >= 5)

if (nrow(context_data) > 0) {
  p_ctx <- ggplot(context_data,
                  aes(x = ctx_up, y = ctx_down, fill = mean_mod)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = n_sites), size = 2.5, colour = "grey20") +
    facet_wrap(~ biotype, nrow = 2) +
    scale_fill_viridis_c(name = "Mean\nmod rate",
                         labels = percent_format(accuracy = 1),
                         option = "plasma") +
    labs(
      title    = "Trinucleotide context of modified cytosines",
      subtitle = "TET+PB condition | numbers = site count | NpCpN context",
      x        = "5' nucleotide (upstream of C)",
      y        = "3' nucleotide (downstream of C)"
    ) +
    theme_srnataps() +
    theme(panel.grid = element_blank())

  save_figure(p_ctx, file.path(opt$figdir, "03e_trinucleotide_context.pdf"),
              width = 12, height = 6)
  message("Figure 3e: Trinucleotide context — done")
}

message("03_modification.R complete.")
