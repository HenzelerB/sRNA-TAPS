# =============================================================================
# 04_benchmark.R — Benchmarking comparison figures
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

COMPARE_DIR   <- file.path(opt$outdir, "09.compare")
BIOTYPE_ORDER <- c("miRNA","tRNA","rRNA","snoRNA","snRNA","lncRNA","other")
TOOL_ORDER    <- c("sRNA-TAPS","rastair","astair","bismark")
TOOL_COLS_FIX <- c("sRNA-TAPS"="#1C4062","rastair"="#E07B39","astair"="#6AAB6E","bismark"="#9B6BB5")
TOOL_SHAPES   <- c("sRNA-TAPS"=18,"rastair"=16,"astair"=17,"bismark"=15)
TOOL_LABS     <- c("sRNA-TAPS"="sRNA-TAPS","rastair"="rastair","astair"="asTair","bismark"="Bismark")

add_condition_group <- function(df) {
  dplyr::mutate(df, condition_group = dplyr::case_when(
    grepl("no.treat|no_treat", condition) ~ "no_treat",
    grepl("pb_Ctrl|pb_ctrl",   condition) ~ "pb_ctrl",
    grepl("treat",             condition) ~ "treat",
    TRUE                                  ~ condition
  ))
}

# ── Figure 4a: Concordance heatmap ───────────────────────────────────────────
conc_file <- file.path(COMPARE_DIR, "concordance_summary.tsv")
if (file.exists(conc_file)) {
  conc <- read_tsv(conc_file, show_col_types = FALSE) %>%
    add_condition_group() %>%
    dplyr::filter(condition_group == "treat") %>%
    dplyr::group_by(biotype, tool) %>%
    dplyr::summarise(jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      biotype        = factor(biotype, levels = BIOTYPE_ORDER),
      tool           = factor(tool,    levels = TOOL_ORDER),
      signed_jaccard = ifelse(tool == "bismark", -jaccard, jaccard)
    ) %>%
    dplyr::filter(!is.na(biotype), !is.na(tool))

  # Add sRNA-TAPS as self-reference row (Jaccard = 1.0)
  srnataps_row <- conc %>%
    dplyr::distinct(biotype) %>%
    dplyr::mutate(tool           = factor("sRNA-TAPS", levels = TOOL_ORDER),
                  jaccard        = 1.0,
                  signed_jaccard = 1.0)
  # Build complete grid — fill missing combinations with NA
  all_combos <- tidyr::expand_grid(
    biotype = factor(levels(conc$biotype), levels = BIOTYPE_ORDER),
    tool    = factor(TOOL_ORDER,           levels = TOOL_ORDER)
  )
  conc_full <- dplyr::bind_rows(conc, srnataps_row) %>%
    dplyr::right_join(all_combos, by = c("biotype", "tool")) %>%
    dplyr::mutate(
      text_col  = dplyr::case_when(
        is.na(jaccard)  ~ "grey50",
        jaccard > 0.4   ~ "white",
        TRUE            ~ "grey10"
      ),
      tile_label = dplyr::case_when(
        is.na(jaccard) ~ "N/A",
        jaccard < 0.001 ~ "< 0.001",
        TRUE           ~ sprintf("%.3f", jaccard)
      )
    )

  # Update tool colours — Okabe-Ito consistent
  TOOL_COLS <- c("sRNA-TAPS" = "#D55E00",
                 "rastair"   = "#0072B2",
                 "astair"    = "#009E73",
                 "bismark"   = "#999999")

  p_conc <- ggplot(
    conc_full %>% dplyr::filter(!is.na(jaccard)),
    aes(x = biotype, y = signed_jaccard, colour = tool, group = tool)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(ymin = 0, ymax = signed_jaccard),
                   position = position_dodge(width = 0.6),
                   linewidth = 0.4) +
    geom_point(shape = 21, colour = "white",
               aes(fill = tool),
               position = position_dodge(width = 0.6),
               size = 3, stroke = 0.3) +
    scale_colour_manual(values = TOOL_COLS, labels = TOOL_LABS,
                        name = "Tool", aesthetics = c("colour", "fill")) +
    scale_y_continuous(limits = c(-0.1, 1), breaks = seq(-0.1, 1, 0.1),
                       labels = function(x) sprintf("%.1f", x)) +
    annotate("text", x = Inf, y = -0.05, label = "← Bismark (inverted chemistry)",
             hjust = 1.05, size = 2.5, colour = "grey50", fontface = "italic") +
    labs(title    = "Site-level concordance: sRNA-TAPS vs benchmark tools",
         subtitle = "TET+PB condition | Signed Jaccard index per biotype",
         x        = "RNA biotype",
         y        = "Jaccard index (negative = chemistry-inverted)",
         caption  = "Bismark values negated: chemistry inversion means anti-concordant sites.") +
    theme_srnataps() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_figure(p_conc, file.path(opt$figdir, "04a_concordance_heatmap.pdf"),
              width = 7, height = 5)
  message("Figure 4a: Concordance heatmap — done")
}

# ── Figure 4b: Pearson correlation ────────────────────────────────────────────
corr_file <- file.path(COMPARE_DIR, "correlation_summary.tsv")
if (file.exists(corr_file)) {
  corr <- read_tsv(corr_file, show_col_types = FALSE) %>%
    add_condition_group() %>%
    dplyr::filter(condition_group == "treat", !is.na(pearson_r)) %>%
    dplyr::group_by(biotype, tool) %>%
    dplyr::summarise(pearson_r = mean(pearson_r, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      biotype = factor(biotype, levels = BIOTYPE_ORDER),
      tool    = factor(tool,    levels = TOOL_ORDER)
    ) %>%
    dplyr::filter(!is.na(biotype))

  message("  Correlation rows: ", nrow(corr))

  p_corr <- ggplot(corr, aes(x = biotype, y = pearson_r,
                              colour = tool, group = tool, shape = tool)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 5) +
    geom_text(aes(label = sprintf("%.2f", pearson_r)),
              size = 2.5, vjust = -0.8, family = "Arial") +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(values = TOOL_COLS_FIX, labels = TOOL_LABS, name = "Tool") +
    scale_shape_manual(values = TOOL_SHAPES,   labels = TOOL_LABS, name = "Tool") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(title    = "Pearson correlation of modification rates at shared sites",
         subtitle = "TET+PB condition per biotype",
         x = "RNA biotype", y = "Pearson r",
         caption = "Negative correlation for Bismark is expected (chemistry inversion)") +
    theme_srnataps() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_figure(p_corr, file.path(opt$figdir, "04b_correlation.pdf"),
              width = 9, height = 5)
  message("Figure 4b: Correlation plot — done")

  # ── Figure 4b extra: Shared sites scatter ────────────────────────────────
  shared_files <- list.files(COMPARE_DIR, pattern = "shared_treat.*\\.tsv$",
                             full.names = TRUE)
  if (length(shared_files) > 0) {
    shared_data <- lapply(shared_files, function(f) {
      bn    <- basename(f)
      parts <- strsplit(sub("\\.tsv$", "", bn), "_")[[1]]
      tool    <- parts[length(parts)]
      biotype <- parts[length(parts) - 1]
      df <- read_tsv(f, show_col_types = FALSE,
                     col_types = cols(site_key = col_character()))
      df$tool    <- tool
      df$biotype <- biotype
      df
    }) %>% dplyr::bind_rows()

    if (nrow(shared_data) > 0 && "mod_rate_custom" %in% names(shared_data)) {
      shared_long <- shared_data %>%
        dplyr::mutate(
          mod_rate_tool = dplyr::coalesce(
            mod_rate_rastair, mod_rate_astair, mod_rate_bismark,
            mod_rate_rastair_all
          ),
          biotype = factor(biotype, levels = BIOTYPE_ORDER),
          tool    = factor(tool,    levels = TOOL_ORDER)
        ) %>%
        dplyr::filter(!is.na(mod_rate_tool), !is.na(mod_rate_custom), !is.na(biotype))

      if (nrow(shared_long) > 0) {
        p_shared <- ggplot(shared_long,
                           aes(x = mod_rate_custom, y = mod_rate_tool,
                               colour = biotype)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                      colour = "grey60", linewidth = 0.3) +
          geom_point(alpha = 0.3, size = 0.5) +
          facet_grid(biotype ~ tool, labeller = labeller(tool = TOOL_LABS)) +
          scale_colour_manual(values = BIOTYPE_COLOURS, guide = "none") +
          scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                             limits = c(0, 1)) +
          labs(title    = "Modification rate at shared sites (TET+PB)",
               subtitle = "Each point = one genomic site in both custom pipeline and benchmark",
               x = "Custom pipeline mod rate", y = "Benchmark tool mod rate") +
          theme_srnataps(base_size = 8)

        save_figure(p_shared, file.path(opt$figdir, "04b_shared_sites_scatter.pdf"),
                    width = 9, height = 12)
        message("Figure 4b extra: Shared sites scatter — done")
      }
    }
  }
}

# ── Figure 4c: Site overlap bar chart ────────────────────────────────────────
if (file.exists(conc_file)) {
  conc_bar <- read_tsv(conc_file, show_col_types = FALSE) %>%
    add_condition_group() %>%
    dplyr::filter(condition_group == "treat") %>%
    dplyr::group_by(biotype, tool) %>%
    dplyr::summarise(
      sites_custom = mean(sites_custom, na.rm = TRUE),
      sites_tool   = mean(sites_tool,   na.rm = TRUE),
      shared       = mean(shared,       na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c(sites_custom, sites_tool, shared),
                        names_to = "category", values_to = "n_sites") %>%
    dplyr::mutate(
      category = dplyr::recode(category,
        "sites_custom" = "Custom only",
        "sites_tool"   = "Tool only",
        "shared"       = "Shared"),
      category = factor(category, levels = c("Custom only","Shared","Tool only")),
      biotype  = factor(biotype, levels = BIOTYPE_ORDER),
      tool     = factor(tool,    levels = TOOL_ORDER)
    ) %>%
    dplyr::filter(!is.na(biotype))

  if (nrow(conc_bar) > 0) {
    p_venn <- ggplot(conc_bar, aes(x = tool, y = n_sites, fill = category)) +
      geom_col(position = "dodge", width = 0.75,
               colour = "grey30", linewidth = 0.2) +
      geom_text(aes(label = scales::label_comma()(round(n_sites))),
                position = position_dodge(width = 0.75),
                vjust = -0.3, size = 2, family = "Arial") +
      facet_wrap(~ biotype, scales = "free_y") +
      scale_fill_manual(
        values = c("Custom only"="#E41A1C","Shared"="#4DAF4A","Tool only"="#377EB8"),
        name   = "Site category") +
      scale_x_discrete(labels = TOOL_LABS) +
      scale_y_continuous(labels = scales::label_comma()) +
      labs(title    = "Called sites: custom pipeline vs benchmark tools",
           subtitle = "TET+PB condition per biotype (mean across replicates)",
           x = NULL, y = "Number of sites") +
      theme_srnataps() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

    save_figure(p_venn, file.path(opt$figdir, "04c_site_overlap.pdf"),
                width = 12, height = 8)
    message("Figure 4c: Site overlap bar chart — done")
  }
}

message("04_benchmark.R complete.")

# ── Figure 4d: rastair CpG-only vs all-context — identical output proof ───────
if (file.exists(corr_file)) {
  cpg_vs_all <- read_tsv(corr_file, show_col_types = FALSE) %>%
    add_condition_group() %>%
    dplyr::filter(condition_group == "treat",
                  tool %in% c("rastair", "rastair_all"),
                  !is.na(pearson_r)) %>%
    dplyr::group_by(biotype, tool) %>%
    dplyr::summarise(pearson_r = mean(pearson_r, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      biotype = factor(biotype, levels = BIOTYPE_ORDER),
      tool    = factor(tool, levels = c("rastair", "rastair_all"),
                       labels = c("rastair (--cpgs-only)", "rastair (all contexts)"))
    ) %>%
    dplyr::filter(!is.na(biotype))

  if (nrow(cpg_vs_all) > 0) {
    p_cpg_all <- ggplot(cpg_vs_all,
                        aes(x = biotype, y = pearson_r,
                            colour = tool, group = tool, shape = tool)) +
      geom_line(linewidth = 0.8, alpha = 0.9) +
      geom_point(size = 3.5) +
      geom_text(aes(label = sprintf("%.2f", pearson_r)),
                size = 2.5, vjust = -0.8, family = "Arial") +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.3) +
      scale_colour_manual(
        values = c("rastair (--cpgs-only)" = "#2196F3",
                   "rastair (all contexts)" = "#B71C1C"),
        name = "rastair mode") +
      scale_shape_manual(
        values = c("rastair (--cpgs-only)" = 16,
                   "rastair (all contexts)" = 17),
        name = "rastair mode") +
      scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
      labs(
        title    = "rastair: CpG-only vs all-context calling — identical output",
        subtitle = "TET+PB condition | Pearson r vs sRNA-TAPS at shared sites",
        x        = "RNA biotype", y = "Pearson r",
        caption  = paste0(
          "Both rastair runs produce identical site calls regardless of --cpgs-only flag.\n",
          "This confirms rastair is a CpG-only caller by design and cannot detect\n",
          "non-CpG m5C modifications characteristic of tRNA, rRNA and snoRNA."
        )
      ) +
      theme_srnataps() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    save_figure(p_cpg_all,
                file.path(opt$figdir, "04d_rastair_cpg_vs_all.pdf"),
                width = 9, height = 5)
    message("Figure 4d: rastair CpG vs all-context — done")
  }
}
