# =============================================================================
# 01_qc.R — QC figures
#
# Inputs (from pipeline output directories):
#   - 02.fastqc/pre_trim/*_fastqc/fastqc_data.txt
#   - 02.fastqc/post_trim/*_fastqc/fastqc_data.txt
#   - logs/align/*_bowtie.log  (or alignment_summary.tsv if generated)
#
# Outputs:
#   - report/figures/01_read_length_distribution.{pdf,png,svg}
#   - report/figures/01_mapping_rates.{pdf,png,svg}
#
# Usage:
#   Rscript 01_qc.R --outdir /path/to/project
# =============================================================================

suppressPackageStartupMessages(library(optparse))
source(file.path(Sys.getenv("SRNATAPS_R_DIR", "/mnt/nfs/home/bhenzeler/projects/RNA_TAPS/sRNA-TAPS/srnataps/report/R"), "00_setup.R"))

option_list <- list(
  make_option("--outdir", type = "character", help = "Project output directory"),
  make_option("--figdir", type = "character", default = NULL,
              help = "Figure output directory (default: outdir/report/figures)")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)


# ── Figure 1a: Read length distribution ──────────────────────────────────────

parse_fastqc_length <- function(fastqc_dir, trim_status) {
  files <- list.files(fastqc_dir, pattern = "fastqc_data.txt",
                      recursive = TRUE, full.names = TRUE)
  all_data <- lapply(files, function(f) {
    sample <- basename(dirname(f))
    sample <- sub("_fastqc$", "", sample)
    sample <- sub("_trimmed$", "", sample)

    lines  <- readLines(f)
    start  <- grep(">>Sequence Length Distribution", lines)
    end    <- grep(">>END_MODULE", lines)

    if (length(start) == 0 || length(end) == 0) return(NULL)
    end  <- end[end > start[1]][1]
    block <- lines[(start[1]+2):(end-1)]

    tryCatch({
      df <- read.table(text = block, sep = "\t", header = FALSE,
                       col.names = c("length_range", "count"))
      df$length    <- as.numeric(sub("-.*", "", df$length_range))
      df$sample    <- sample
      df$condition <- get_condition(sample)
      df$cell_line <- get_cell_line(sample)
      df$trim      <- trim_status
      df
    }, error = function(e) NULL)
  })
  dplyr::bind_rows(Filter(Negate(is.null), all_data))
}

pre_dir  <- file.path(opt$outdir, "02.fastqc", "pre_trim")
post_dir <- file.path(opt$outdir, "02.fastqc", "post_trim")

if (dir.exists(pre_dir) && dir.exists(post_dir)) {
  pre  <- parse_fastqc_length(pre_dir,  "Pre-trim")
  post <- parse_fastqc_length(post_dir, "Post-trim")
  len_data <- dplyr::bind_rows(pre, post) %>%
    dplyr::mutate(
      condition = factor(condition, levels = names(CONDITION_COLOURS)),
      trim      = factor(trim, levels = c("Pre-trim", "Post-trim"))
    )

  p_len <- ggplot(len_data, aes(x = length, y = count,
                                colour = condition, group = sample)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    facet_grid(cell_line ~ trim) +
    scale_colour_manual(values = CONDITION_COLOURS, labels = CONDITION_LABELS,
                        name = "Condition") +
    scale_y_log10(labels = label_comma()) +
    scale_x_continuous(breaks = c(15, 20, 25, 30, 40, 50)) +
    labs(
      title    = "Read length distribution",
      subtitle = "Per sample, pre- and post-adapter trimming",
      x        = "Read length (nt)",
      y        = "Read count (log10)",
      caption  = "TruSeq small RNA adapter trimmed with Trim Galore --small_rna"
    ) +
    theme_srnataps()

  save_figure(p_len, file.path(opt$figdir, "01a_read_length_distribution.pdf"),
              width = 9, height = 6)
  message("Figure 1a: Read length distribution — done")
} else {
  message("WARNING: FastQC directories not found — skipping Figure 1a")
}


# ── Figure 1b: Mapping rates ──────────────────────────────────────────────────

# Parse bowtie logs: extract alignment rate per sample
parse_bowtie_logs <- function(log_dir) {
  files <- list.files(log_dir, pattern = "_bowtie\\.log$", full.names = TRUE)
  if (length(files) == 0) return(NULL)

  lapply(files, function(f) {
    sample <- sub("_bowtie\\.log$", "", basename(f))
    lines  <- readLines(f)

    total_line  <- grep("reads processed", lines, value = TRUE)
    mapped_line <- grep("alignment rate",  lines, value = TRUE)
    multi_line  <- grep("reported alignments", lines, value = TRUE)

    total  <- as.numeric(gsub("[^0-9]", "", total_line[1]))
    rate   <- as.numeric(gsub("[^0-9\\.]", "", mapped_line[1]))

    if (is.na(total) || is.na(rate)) return(NULL)

    data.frame(
      sample       = sample,
      total_reads  = total,
      mapping_rate = rate,
      condition    = get_condition(sample),
      cell_line    = get_cell_line(sample),
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
}

log_dir  <- file.path(opt$outdir, "logs", "align")
map_data <- parse_bowtie_logs(log_dir)

if (!is.null(map_data) && nrow(map_data) > 0) {
  map_data <- map_data %>%
    dplyr::arrange(condition, cell_line, sample) %>%
    dplyr::mutate(
      sample    = factor(sample, levels = unique(sample)),
      condition = factor(condition, levels = names(CONDITION_COLOURS))
    )

  p_map <- ggplot(map_data, aes(x = sample, y = mapping_rate, fill = condition)) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = 60, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    facet_wrap(~ cell_line, scales = "free_x") +
    scale_fill_manual(values = CONDITION_COLOURS, labels = CONDITION_LABELS,
                      name = "Condition") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Bowtie1 alignment rates",
      subtitle = "Percentage of trimmed reads aligning to hg38",
      x        = NULL,
      y        = "Alignment rate (%)",
      caption  = "Dashed line: 60% reference threshold. Parameters: -v2 --norc -k10 --best --strata -m100"
    ) +
    theme_srnataps() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  save_figure(p_map, file.path(opt$figdir, "01b_mapping_rates.pdf"),
              width = 10, height = 5)
  message("Figure 1b: Mapping rates — done")
} else {
  message("WARNING: Bowtie logs not found — skipping Figure 1b")
}

message("01_qc.R complete.")
