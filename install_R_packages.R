#!/usr/bin/env Rscript
# =============================================================================
# install_R_packages.R — Install all R packages needed for sRNA-TAPS reports
#
# Usage:
#   Rscript install_R_packages.R
#
# Run once after setting up the conda environment.
# Installs to the default R library path.
# =============================================================================

pkgs <- c(
  "ggplot2",      # plotting
  "dplyr",        # data manipulation
  "tidyr",        # reshaping
  "readr",        # fast file reading
  "scales",       # axis formatting
  "RColorBrewer", # colour palettes
  "ggrepel",      # non-overlapping labels
  "patchwork",    # combining plots
  "cowplot",      # publication themes
  "viridis",      # colour scales
  "ggbeeswarm",   # beeswarm plots
  "ggridges",     # ridge plots
  "optparse",     # CLI argument parsing
  "stringr"       # string manipulation
)

installed <- rownames(installed.packages())
to_install <- pkgs[!pkgs %in% installed]

if (length(to_install) > 0) {
  message("Installing ", length(to_install), " packages: ",
          paste(to_install, collapse = ", "))
  install.packages(to_install,
                   repos   = "https://cloud.r-project.org",
                   quiet   = TRUE,
                   dependencies = TRUE)
} else {
  message("All R packages already installed.")
}

# Verify
ok      <- c()
missing <- c()
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    ok <- c(ok, p)
  } else {
    missing <- c(missing, p)
  }
}

message("\n=== R package status ===")
for (p in ok)      message("  OK      ", p)
for (p in missing) message("  MISSING ", p)

if (length(missing) > 0) {
  stop("Some packages failed to install: ", paste(missing, collapse = ", "))
} else {
  message("\nAll R packages ready.")
}

# ── Bioconductor packages (for sequence logos) ────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "Biostrings")
bioc_installed <- rownames(installed.packages())
bioc_to_install <- bioc_pkgs[!bioc_pkgs %in% bioc_installed]

if (length(bioc_to_install) > 0) {
  message("Installing Bioconductor packages: ",
          paste(bioc_to_install, collapse = ", "))
  BiocManager::install(bioc_to_install, ask = FALSE, update = FALSE)
}

# ggseqlogo
if (!"ggseqlogo" %in% rownames(installed.packages())) {
  install.packages("ggseqlogo",
                   repos = "https://cloud.r-project.org",
                   quiet = TRUE)
}

message("Bioconductor packages done.")
