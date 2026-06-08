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

# ── Colour scheme (aligned with sRNA-TAPS logo gradient) ─────────────────────
# Logo palette: dark navy #1C4062 → primary blue #3B7DB6 → teal #3898BD → light #7CBBD4
SRNATAPS_NAVY  <- "#1C4062"
SRNATAPS_BLUE  <- "#3B7DB6"
SRNATAPS_MID   <- "#3683AF"
SRNATAPS_TEAL  <- "#3898BD"
SRNATAPS_LIGHT <- "#7CBBD4"

COLS <- c(
  treat    = "#1C4062",
  pb_Ctrl  = "#3B7DB6",
  no_treat = "#7CBBD4",
  HEK      = "#3898BD",
  Caco2    = "#3683AF",
  high     = "#1C4062",
  medium   = "#3B7DB6",
  low      = "#A8CADA"
)

# Condition colours — keyed to internal condition strings
CONDITION_COLOURS <- c(
  "treat"    = "#1C4062",
  "pb_ctrl"  = "#3B7DB6",
  "no_treat" = "#7CBBD4",
  "old"      = "#C8DCE8"
)

CONDITION_LABELS <- c(
  "treat"    = "TET + PB",
  "pb_ctrl"  = "PB only",
  "no_treat" = "Untreated",
  "old"      = "Old HEK"
)

# Cell line colours
CELL_COLOURS <- c(
  "HEK"   = "#1C4062",
  "Caco2" = "#3898BD"
)

# Cell line shapes
CELL_SHAPES <- c(
  "HEK"   = 16,
  "Caco2" = 17
)

# Confidence colours
CONF_COLOURS <- c(
  "High"   = "#1C4062",
  "Medium" = "#3B7DB6",
  "Low"    = "#A8CADA"
)

# Biotype colours — sequential from logo palette + complementary neutrals
BIOTYPE_COLOURS <- c(
  "miRNA"  = "#1C4062",
  "tRNA"   = "#3B7DB6",
  "rRNA"   = "#3683AF",
  "snoRNA" = "#3898BD",
  "snRNA"  = "#7CBBD4",
  "piRNA"  = "#A8CADA",
  "lncRNA" = "#5B8FA8",
  "other"  = "#8DAFC0"
)

# Benchmark tool colours — sRNA-TAPS gets logo primary blue; others get neutral tones
TOOL_COLOURS_BENCH <- c(
  "sRNA-TAPS" = "#1C4062",
  "rastair"   = "#E07B39",
  "astair"    = "#6AAB6E",
  "bismark"   = "#9B6BB5"
)

# ── Shared Arial 8pt theme ────────────────────────────────────────────────────
theme_arial8 <- theme_bw(base_size = 8, base_family = "Arial") +  # sRNA-TAPS house style: Arial 8pt
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
