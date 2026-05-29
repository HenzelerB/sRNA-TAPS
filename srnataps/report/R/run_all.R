#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — Run all sRNA-TAPS report figure scripts
#
# Usage:
#   Rscript run_all.R --outdir /path/to/project [--figdir /path/to/figures]
#
# Runs in order:
#   01_qc.R           — Read length distribution, mapping rates
#   02_biotype.R      — Biotype composition
#   03_modification.R — Modification rates, top species, conditions
#   04_benchmark.R    — Benchmarking concordance and correlation
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--outdir",   type = "character", help = "Project output directory"),
  make_option("--figdir",   type = "character", default = NULL,
              help = "Figure output directory (default: outdir/report/figures)"),
  make_option("--scripts",  type = "character", default = NULL,
              help = "Directory containing R scripts (default: same dir as run_all.R)"),
  make_option("--skip-qc",  action = "store_true", default = FALSE),
  make_option("--skip-bio", action = "store_true", default = FALSE),
  make_option("--skip-mod", action = "store_true", default = FALSE),
  make_option("--skip-bench", action = "store_true", default = FALSE),
  make_option("--skip-logos", action = "store_true", default = FALSE,
              help = "Skip sequence logo generation (requires BSgenome)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$outdir)) stop("--outdir is required")
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
if (is.null(opt$scripts)) opt$scripts <- dirname(sys.frame(1)$ofile)

dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

run_script <- function(script_name, extra_args = "") {
  script_path <- file.path(opt$scripts, script_name)
  if (!file.exists(script_path)) {
    message("WARNING: Script not found: ", script_path)
    return(invisible(NULL))
  }
  cmd <- paste(
    "Rscript", shQuote(script_path),
    "--outdir", shQuote(opt$outdir),
    "--figdir", shQuote(opt$figdir),
    extra_args
  )
  message("\n", strrep("=", 60))
  message("Running: ", script_name)
  message(strrep("=", 60))
  t_start <- proc.time()
  result  <- system(cmd)
  elapsed <- round((proc.time() - t_start)[["elapsed"]], 1)
  if (result != 0) {
    message("WARNING: ", script_name, " exited with code ", result)
  } else {
    message("Done in ", elapsed, "s")
  }
  invisible(result)
}

message("sRNA-TAPS report generation")
message("Output directory : ", opt$outdir)
message("Figure directory : ", opt$figdir)
message(strrep("-", 60))

if (!opt$`skip-qc`)    run_script("01_qc.R")
if (!opt$`skip-bio`)   run_script("02_biotype.R")
if (!opt$`skip-mod`)   run_script("03_modification.R")
if (!opt$`skip-bench`) run_script("04_benchmark.R")
if (!opt$`skip-logos`) run_script("05_sequence_logos.R")

message("\n", strrep("=", 60))
message("All figures complete.")
message("Output: ", opt$figdir)
message(strrep("=", 60))

# Print summary of files produced
figs <- list.files(opt$figdir, pattern = "\\.(pdf|png|svg)$",
                   full.names = FALSE)
message("\nFigures produced (", length(figs), "):")
for (f in sort(figs)) message("  ", f)
