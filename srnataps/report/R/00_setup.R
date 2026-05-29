# =============================================================================
# 00_setup.R — Shared theme, colours, and helper functions
# Sourced by all other R scripts in the report pipeline.
# Colour scheme and theme from akschneider analysis pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(viridis)
  library(ggbeeswarm)
  library(ggridges)
  library(extrafont)
  library(stringr)
})

# Load Arial font for PDF output
tryCatch(
  loadfonts(device = "pdf", quiet = TRUE),
  error = function(e) message("Note: extrafont loadfonts failed — using sans serif fallback")
)

# ── Colour scheme ─────────────────────────────────────────────────────────────
COLS <- c(
  treat    = "#C0392B",
  pb_Ctrl  = "#E67E22",
  no_treat = "#2980B9",
  HEK      = "#8E44AD",
  Caco2    = "#27AE60",
  high     = "#C0392B",
  medium   = "#E67E22",
  low      = "#7F8C8D"
)

# Condition colours — keyed to internal condition strings
CONDITION_COLOURS <- c(
  "treat"    = "#C0392B",
  "pb_ctrl"  = "#E67E22",
  "no_treat" = "#2980B9",
  "old"      = "#BDC3C7"
)

CONDITION_LABELS <- c(
  "treat"    = "TET + PB",
  "pb_ctrl"  = "PB only",
  "no_treat" = "Untreated",
  "old"      = "Old HEK"
)

# Cell line colours
CELL_COLOURS <- c(
  "HEK"   = "#8E44AD",
  "Caco2" = "#27AE60"
)

# Cell line shapes
CELL_SHAPES <- c(
  "HEK"   = 16,
  "Caco2" = 17
)

# Confidence colours — capitalised to match data values
CONF_COLOURS <- c(
  "High"   = "#C0392B",
  "Medium" = "#E67E22",
  "Low"    = "#7F8C8D"
)

# Biotype colours
BIOTYPE_COLOURS <- c(
  "miRNA"  = "#C0392B",
  "tRNA"   = "#2980B9",
  "rRNA"   = "#27AE60",
  "snoRNA" = "#8E44AD",
  "snRNA"  = "#E67E22",
  "piRNA"  = "#16A085",
  "lncRNA" = "#F39C12",
  "other"  = "#7F8C8D"
)

# ── Shared Arial 8pt theme ────────────────────────────────────────────────────
theme_arial8 <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    text               = element_text(size = 8, family = "Arial", colour = "black"),
    axis.text          = element_text(size = 8, family = "Arial", colour = "black"),
    axis.title         = element_text(size = 8, family = "Arial", colour = "black", face = "bold"),
    plot.title         = element_text(size = 8, family = "Arial", colour = "black"),
    plot.subtitle      = element_text(size = 6, family = "Arial", colour = "black", face = "italic"),
    legend.text        = element_text(size = 8, family = "Arial", colour = "black"),
    legend.title       = element_text(size = 8, family = "Arial", colour = "black"),
    strip.text         = element_text(size = 8, family = "Arial", colour = "black", face = "bold"),
    strip.background   = element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    legend.position    = "right",
    legend.key.size    = unit(0.35, "cm"),
    plot.margin        = margin(4, 4, 4, 4)
  )

# Alias so other scripts can call theme_srnataps() as well
theme_srnataps <- function(...) theme_arial8

# ── Save helper ───────────────────────────────────────────────────────────────
# Always saves PDF (vector, for publication) + PNG (raster) + SVG
save_figure <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  base    <- tools::file_path_sans_ext(filename)
  pdf_out <- paste0(base, ".pdf")
  png_out <- paste0(base, ".png")
  svg_out <- paste0(base, ".svg")

  ggsave(pdf_out, plot = plot, width = width, height = height,
         device = cairo_pdf, family = "Arial")
  ggsave(png_out, plot = plot, width = width, height = height, dpi = dpi)
  ggsave(svg_out, plot = plot, width = width, height = height, device = "svg")

  message("  Saved: ", pdf_out)
}

# ── Sample metadata helpers ───────────────────────────────────────────────────
get_condition <- function(sample) {
  dplyr::case_when(
    grepl("pb_Ctrl",  sample) ~ "pb_ctrl",
    grepl("no-treat", sample) ~ "no_treat",
    grepl("treat",    sample) ~ "treat",
    grepl("old",      sample) ~ "old",
    TRUE                      ~ "unknown"
  )
}

get_cell_line <- function(sample) {
  dplyr::case_when(
    grepl("HEK",   sample) ~ "HEK",
    grepl("Caco2", sample) ~ "Caco2",
    TRUE                   ~ "unknown"
  )
}

get_replicate <- function(sample) {
  stringr::str_extract(sample, "R[0-9]+")
}

message("sRNA-TAPS theme loaded (Arial 8pt).")
