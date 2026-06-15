# =============================================================================
# 00_setup.R — Shared theme, colours, and helper functions
# Sourced by all other R scripts in the report pipeline.
#
# GENERIC DESIGN: reads samples.tsv from --outdir to build sample metadata
# lookups dynamically. Works for any cell lines, any condition names, any
# number of replicates — not just HEK/Caco2 or TAPS-specific condition names.
#
# Falls back gracefully to sample-name string matching when samples.tsv is not
# found, preserving backward compatibility.
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


# =============================================================================
# 1. Locate samples.tsv via --outdir arg passed to the parent script
# =============================================================================

.get_outdir_from_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- which(args == "--outdir")
  if (length(idx) > 0 && idx[1] < length(args)) return(args[idx[1] + 1])
  env <- Sys.getenv("SRNATAPS_OUTDIR", unset = "")
  if (nchar(env) > 0) return(env)
  getwd()
}

SRNATAPS_OUTDIR   <- .get_outdir_from_args()
SRNATAPS_SAMPLES  <- file.path(SRNATAPS_OUTDIR, "samples.tsv")


# =============================================================================
# 2. Read samples.tsv and build lookup tables
# =============================================================================

if (file.exists(SRNATAPS_SAMPLES)) {
  .meta <- readr::read_tsv(SRNATAPS_SAMPLES, show_col_types = FALSE)

  # Normalise condition strings: replace hyphens/spaces with underscores
  .meta <- .meta %>%
    dplyr::mutate(
      condition_norm = gsub("[-. ]", "_", condition),
      cell_line      = as.character(cell_line)
    )

  # Named vectors for fast O(1) lookup: sample → condition / cell_line
  .COND_LOOKUP  <- setNames(.meta$condition_norm, .meta$sample)
  .CL_LOOKUP    <- setNames(.meta$cell_line,      .meta$sample)

  CONDITIONS  <- unique(.meta$condition_norm)
  CELL_LINES  <- unique(.meta$cell_line)
  message("  Loaded ", nrow(.meta), " samples | ",
          length(CONDITIONS), " conditions | ",
          length(CELL_LINES), " cell lines")
} else {
  message("  WARNING: samples.tsv not found at ", SRNATAPS_SAMPLES,
          " — falling back to sample-name parsing")
  .COND_LOOKUP <- NULL
  .CL_LOOKUP   <- NULL
  CONDITIONS   <- c("treat", "pb_ctrl", "no_treat")
  CELL_LINES   <- c("HEK", "Caco2")
}


# =============================================================================
# 3. Helper functions — use lookup table when available, name-parse as fallback
# =============================================================================

get_condition <- function(sample) {
  if (!is.null(.COND_LOOKUP)) {
    res <- .COND_LOOKUP[sample]
    res[is.na(res)] <- "unknown"
    return(unname(res))
  }
  # Fallback: infer from sample name
  dplyr::case_when(
    grepl("pb_Ctrl|pb_ctrl", sample) ~ "pb_ctrl",
    grepl("no[-_]treat",     sample) ~ "no_treat",
    grepl("treat",           sample) ~ "treat",
    grepl("old",             sample) ~ "old",
    TRUE                             ~ "unknown"
  )
}

get_cell_line <- function(sample) {
  if (!is.null(.CL_LOOKUP)) {
    res <- .CL_LOOKUP[sample]
    res[is.na(res)] <- "unknown"
    return(unname(res))
  }
  # Fallback: infer from sample name
  dplyr::case_when(
    grepl("HEK",   sample) ~ "HEK",
    grepl("Caco2", sample) ~ "Caco2",
    TRUE                   ~ "unknown"
  )
}

get_replicate <- function(sample) {
  stringr::str_extract(sample, "R[0-9]+")
}


# =============================================================================
# 4. sRNA-TAPS colour palette (logo gradient: navy → teal)
# =============================================================================

SRNATAPS_NAVY  <- "#1C4062"
SRNATAPS_BLUE  <- "#3B7DB6"
SRNATAPS_MID   <- "#3683AF"
SRNATAPS_TEAL  <- "#3898BD"
SRNATAPS_LIGHT <- "#7CBBD4"
SRNATAPS_PALE  <- "#A8CADA"

# Base palette in order (6 steps from dark to light)
.BASE_PALETTE <- c(SRNATAPS_NAVY, SRNATAPS_BLUE, SRNATAPS_MID,
                   SRNATAPS_TEAL, SRNATAPS_LIGHT, SRNATAPS_PALE)

# Generate n colours from the sRNA-TAPS gradient
srnataps_palette <- function(n) {
  if (n == 0) return(character(0))
  if (n <= length(.BASE_PALETTE)) return(.BASE_PALETTE[seq_len(n)])
  colorRampPalette(.BASE_PALETTE)(n)
}


# =============================================================================
# 5. Dynamic condition colours and labels
# =============================================================================

# Known display labels for common TAPS condition names.
# Add your own by setting SRNATAPS_CONDITION_LABELS env var as
# "cond1=Label 1,cond2=Label 2"
.KNOWN_LABELS <- c(
  "treat"    = "PB⁺ & TET⁺",
  "pb_ctrl"  = "PB⁺",
  "no_treat" = "PB⁻ & TET⁻",
  "no-treat" = "PB⁻ & TET⁻"
)

# Parse any user-supplied overrides from environment
.env_labels <- Sys.getenv("SRNATAPS_CONDITION_LABELS", unset = "")
if (nchar(.env_labels) > 0) {
  .pairs <- strsplit(strsplit(.env_labels, ",")[[1]], "=")
  .user  <- setNames(
    sapply(.pairs, `[`, 2),
    sapply(.pairs, `[`, 1)
  )
  .KNOWN_LABELS <- c(.user, .KNOWN_LABELS[!names(.KNOWN_LABELS) %in% names(.user)])
}

# Build condition colours: palette steps assigned in the order conditions appear
# ── CANONICAL FIGURE PALETTE (Aurora Borealis, brightened) ───────────────────
# Single source of truth for all report figures. Condition colours below;
# biotype spectrum (BIOTYPE_COLOURS) is interpolated from the same three anchors.
# Known TAPS conditions get fixed colours; custom sets fall back to srnataps_palette.
.KNOWN_CONDITION_COLOURS <- c(
  "treat"    = "#FFD680",   # light gold  — TET+PB
  "pb_ctrl"  = "#5BAFD0",   # light blue  — PB
  "no_treat" = "#4DDEB8"    # light mint  — UNTREATED
)
CONDITION_COLOURS <- setNames(
  sapply(CONDITIONS, function(c) {
    if (c %in% names(.KNOWN_CONDITION_COLOURS)) .KNOWN_CONDITION_COLOURS[[c]]
    else srnataps_palette(length(CONDITIONS))[which(CONDITIONS == c)]
  }),
  CONDITIONS
)

# Build condition labels: use known label if available, else use condition name
CONDITION_LABELS <- setNames(
  sapply(CONDITIONS, function(c) {
    if (c %in% names(.KNOWN_LABELS)) .KNOWN_LABELS[[c]] else c
  }),
  CONDITIONS
)


# =============================================================================
# 6. Dynamic cell line colours and shapes
# =============================================================================

.CL_SHAPES_BASE <- c(16, 17, 15, 18, 3, 4, 8, 1)   # filled circle, triangle, ...

CELL_COLOURS <- setNames(
  srnataps_palette(length(CELL_LINES)),
  CELL_LINES
)

CELL_SHAPES <- setNames(
  .CL_SHAPES_BASE[seq_along(CELL_LINES)],
  CELL_LINES
)


# =============================================================================
# 7. Fixed palettes (biotypes + benchmark tools — pipeline-defined, not dynamic)
# =============================================================================

# Okabe-Ito colorblind-safe palette (8 colours)
# Biotype colours: spectrum derived from Aurora Borealis condition palette
BIOTYPE_COLOURS <- setNames(
  colorRampPalette(c("#FFD680", "#5BAFD0", "#4DDEB8"))(8),
  c("miRNA", "tRNA", "rRNA", "snoRNA", "snRNA", "piRNA", "lncRNA", "other")
)

TOOL_COLOURS_BENCH <- c(
  "sRNA-TAPS" = "#1C4062",
  "rastair"   = "#E07B39",
  "astair"    = "#6AAB6E",
  "bismark"   = "#9B6BB5"
)

CONF_COLOURS <- c(
  "High"   = "#1C4062",
  "Medium" = "#3B7DB6",
  "Low"    = "#A8CADA"
)


# =============================================================================
# 8. Theme (Arial 8pt — sRNA-TAPS house style)
# =============================================================================

theme_arial8 <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    text               = element_text(size = 8, family = "Arial", colour = "black"),
    axis.text          = element_text(size = 8, family = "Arial", colour = "black"),
    axis.title         = element_text(size = 8, family = "Arial", colour = "black",
                                      face = "bold"),
    plot.title         = element_text(size = 8, family = "Arial", colour = "black"),
    plot.subtitle      = element_text(size = 6, family = "Arial", colour = "grey50",
                                      face = "plain"),
    legend.text        = element_text(size = 8, family = "Arial", colour = "black"),
    legend.title       = element_text(size = 8, family = "Arial", colour = "black"),
    strip.text         = element_text(size = 8, family = "Arial", colour = "black",
                                      face = "bold"),
    strip.background   = element_rect(fill = "grey95", colour = "grey60",
                                      linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    legend.position    = "right",
    legend.key.size    = unit(0.35, "cm"),
    plot.margin        = margin(4, 4, 4, 4)
  )

# Alias
theme_srnataps <- function(...) theme_arial8


# =============================================================================
# 9. Save helper — PDF + PNG + SVG
# =============================================================================

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


# =============================================================================
# 10. Summary
# =============================================================================

message("sRNA-TAPS theme loaded (Arial 8pt).")
message("  Conditions : ", paste(CONDITIONS,  collapse = ", "))
message("  Cell lines : ", paste(CELL_LINES,  collapse = ", "))
# ── Auto-detect TAPS condition roles ─────────────────────────────────────────
# Matches standard TAPS naming (treat / pb* / no_treat).
# Override via env vars: SRNATAPS_TREAT_COND, SRNATAPS_CTRL_COND, SRNATAPS_UNTR_COND
.ev_treat <- Sys.getenv("SRNATAPS_TREAT_COND", unset = "")
.ev_ctrl  <- Sys.getenv("SRNATAPS_CTRL_COND",  unset = "")
.ev_untr  <- Sys.getenv("SRNATAPS_UNTR_COND",  unset = "")

TREAT_COND <- if (nchar(.ev_treat) > 0) .ev_treat else {
  c1 <- CONDITIONS[CONDITIONS == "treat"][1]
  if (!is.na(c1)) c1 else
    CONDITIONS[grepl("treat", CONDITIONS) & !grepl("^no", CONDITIONS)][1]
}
CTRL_COND <- if (nchar(.ev_ctrl) > 0) .ev_ctrl else {
  c1 <- CONDITIONS[grepl("^pb", CONDITIONS) & !grepl("treat", CONDITIONS)][1]
  if (!is.na(c1)) c1 else
    CONDITIONS[grepl("ctrl|control", CONDITIONS, ignore.case=TRUE) & !grepl("treat", CONDITIONS)][1]
}
UNTR_COND <- if (nchar(.ev_untr) > 0) .ev_untr else {
  c1 <- CONDITIONS[grepl("^no_?treat|^untr", CONDITIONS)][1]
  if (!is.na(c1)) c1 else CONDITIONS[!CONDITIONS %in% c(TREAT_COND, CTRL_COND)][1]
}

# Final fallback
if (is.na(TREAT_COND)) TREAT_COND <- CONDITIONS[length(CONDITIONS)]
if (is.na(CTRL_COND))  CTRL_COND  <- CONDITIONS[min(2, length(CONDITIONS))]
if (is.na(UNTR_COND))  UNTR_COND  <- CONDITIONS[1]

message("  Treat cond : ", TREAT_COND)
message("  Ctrl cond  : ", CTRL_COND)
message("  Untr cond  : ", UNTR_COND)

