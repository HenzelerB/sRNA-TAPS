# =============================================================================
# 05_sequence_logos.R — sequence context around modified cytosines
#
# Changes in this version (2026-10-03):
#
#   --site-table <tsv> (default unset = no site filter):
#   keep only foreground and background sites with strict_canonical == True
#   in the shared site-classification table, for EVERY biotype. Rows are
#   matched on (biotype_folder, cell_line, chrom, start, strand) <->
#   (biotype, cell_line, chrom, start, strand), per cell line, before
#   replicate collapsing and before cell-line pooling. Sites absent from the
#   table are dropped (counted in the log and in
#   <figdir>/06_sequence_context_site_table_match.tsv).
#
# Two further fixes:
#
#   A. Chromosome naming. normalise_chrom() used to turn Ensembl contig names
#      such as KI270733.1 into "chrKI270733.1", which is not a BSgenome
#      seqname (hg38: chrUn_KI270733v1), so extract_sequences() silently
#      dropped every site on unplaced/unlocalized contigs (incl. the 45S rDNA
#      contigs GL000220.1, KI270733.1, KI270442.1). Contig names are now
#      resolved by lookup in seqnames(genome) (chrUn_*, chrN_*_random, ...).
#      Unresolvable names and getSeq() failures are errors; windows dropped at
#      contig ends or for containing N are counted in the log.
#   B. Test. Per-position p-values now come from chisq.test(...,
#      simulate.p.value = TRUE, B = 100000) with set.seed(42) once per panel
#      before the position loop (uniformly for all biotypes), then BH across
#      the window as before. The chi-squared statistic and min_expected are
#      still reported. All-zero base columns are dropped before testing.
#
# Changes in the 2026-10-02 version:
#
#   1. Double-bind bug fixed. The original's
#        lst %>% dplyr::bind_rows(Filter(Negate(is.null), .))
#      expands to bind_rows(lst, Filter(Negate(is.null), lst)) (magrittr also
#      passes `.` as the first argument), so every foreground and background
#      row was included twice. Lists are now assigned first, then bound once.
#   2. Replicate rows collapsed to unique genomic sites, so each C counts once:
#        Foreground: a site (cell_line, biotype, chrom, start, strand) is kept
#          once if it passes snp_flag == PASS, coverage >= min_cov,
#          mod_rate >= min_mod and padj < max_padj in at least one replicate.
#        Background: unique sites passing snp_flag == PASS,
#          coverage >= min_cov and mod_rate <= 0.05 in at least one replicate,
#          EXCLUDING any site that is in the foreground for the same
#          cell_line and biotype.
#        Background is then sampled (set.seed(42)) per panel to
#          max(3 x foreground n, 500) unique sites (the original's 3x rule and
#          500 floor), or all available sites if fewer exist.
#   3. --hw-c (default 5): half-window of the enrichment panel (06c),
#      previously hardcoded as hw_c <- 5.
#   4. The per-position chi-squared table is always written to
#      <figdir>/06c_sequence_context_enrichment_significance.tsv, and the
#      per-base foreground/background frequencies to
#      <figdir>/06c_sequence_context_enrichment_frequencies.tsv.
#   5. --pool-cell-lines (default off): one 06c panel per biotype, with Caco2
#      and HEK unique sites unioned (a site modified in both lines counts
#      once; background excludes any site in the pooled foreground).
#      Default (off) = one panel per biotype x cell_line.
#   6. Layout: shorter panel titles aligned to the whole plot (no clipping),
#      shorter y title ("Δ freq. (pp)") that fits the panel height,
#      fixed panels of 70 x 18 mm at +/-5 nt and 85 x 18 mm above (original:
#      55 x 14.9 mm), every other x tick
#      labelled when hw_c > 5, and a two-column layout (Caco2 left, HEK
#      right) for per-cell-line output when hw_c > 5.
#   7. Background criterion is mod_rate <= 0.05 (original: < 0.05), matching
#      the manuscript's stated definition.
#
# Requires Bioconductor packages:
#   BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges"))
#   install.packages(c("ggseqlogo", "gridExtra"))
#
# Inputs:
#   - 07.taps_calls/<biotype>/<sample>_<biotype>_taps.tsv
#     Required columns: chrom, start, strand, coverage, mod_rate, padj, snp_flag
#
# Outputs (into --figdir):
#   - 06a_sequence_context_all.{pdf,png,svg}          (fixed +/-10 nt logo)
#   - 06b_sequence_context_by_biotype.{pdf,png,svg}   (fixed +/-10 nt logos)
#   - 06c_sequence_context_enrichment.{pdf,png,svg}   (+/- --hw-c nt)
#   - 06c_sequence_context_enrichment_significance.tsv
#   - 06c_sequence_context_enrichment_frequencies.tsv
#   - 06_sequence_context_site_summary.tsv            (unique fg/bg n per panel)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggseqlogo)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(stringr)
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
  make_option("--outdir",    type = "character", help = "Project output directory"),
  make_option("--figdir",    type = "character", default = NULL),
  # These are supplied by the workflow.  --genome is retained for interface
  # compatibility; sequence extraction continues to use the installed hg38
  # BSgenome object below.
  make_option("--calls-dir", type = "character", default = NULL,
              help = "TAPS calls directory (defaults to <outdir>/07.taps_calls)"),
  make_option("--genome",    type = "character", default = NULL,
              help = "Reference FASTA (accepted for workflow compatibility)"),
  make_option("--biotypes",  type = "character",
              default = "miRNA,tRNA,rRNA,snoRNA",
              help = "Comma-separated biotypes [default: miRNA,tRNA,rRNA,snoRNA]"),
  make_option("--min-cov",   type = "integer",   default = 5),
  make_option("--min-mod",   type = "double",    default = 0.50,
              help = "Minimum mod_rate for high-confidence sites [default: 0.50]"),
  make_option("--max-padj",  type = "double",    default = 0.01,
              help = "Maximum adjusted p-value [default: 0.01]"),
  make_option("--max-sites", type = "integer",   default = 2000,
              help = "Max sites per 06a/06b logo (random sample if more) [default: 2000]"),
  make_option("--condition", type = "character", default = "treat",
              help = "Condition to use [default: treat]"),
  make_option("--hw-c",      type = "integer",   default = 5,
              help = "Half-window (nt) of the background-corrected enrichment panel (06c) [default: 5]"),
  make_option("--bg-max-mod", type = "double",   default = 0.05,
              help = "Maximum mod_rate for background (unmodified) C [default: 0.05]"),
  make_option("--pool-cell-lines", action = "store_true", default = FALSE,
              help = "One 06c panel per biotype with cell lines unioned [default: per cell line]"),
  make_option("--site-table", type = "character", default = NULL,
              help = paste0("Site-classification TSV; keep only sites with strict_canonical == True, ",
                            "joined on (biotype_folder, cell_line, chrom, start, strand) [default: unset]"))
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$figdir)) opt$figdir <- file.path(opt$outdir, "report", "figures")
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

BIOTYPES_USE <- strsplit(opt$biotypes, ",")[[1]]
CALLS_DIR    <- if (!is.null(opt$`calls-dir`)) {
  opt$`calls-dir`
} else {
  file.path(opt$outdir, "07.taps_calls")
}
if (!dir.exists(CALLS_DIR)) {
  legacy_calls_dir <- file.path(opt$outdir, "07.taps_calls_3letter")
  if (dir.exists(legacy_calls_dir)) CALLS_DIR <- legacy_calls_dir
}
GENOME       <- BSgenome.Hsapiens.UCSC.hg38
WINDOWS      <- 10
POOL         <- isTRUE(opt$`pool-cell-lines`)
POOLED_LABEL <- "Caco2 + HEK"
SITE_KEYS    <- c("biotype", "cell_line", "chrom", "start", "strand")

message("Options: hw-c = ", opt$`hw-c`, " | pool-cell-lines = ", POOL,
        " | bg mod_rate <= ", opt$`bg-max-mod`,
        " | site-table = ", if (is.null(opt$`site-table`)) "unset" else opt$`site-table`)

# ── --site-table: strict canonical site set ──────────────────────────────────
SITE_TABLE <- NULL
if (!is.null(opt$`site-table`)) {
  if (!file.exists(opt$`site-table`)) stop("--site-table not found: ", opt$`site-table`)
  SITE_TABLE <- readr::read_tsv(
    opt$`site-table`, show_col_types = FALSE,
    col_select = c("biotype_folder", "cell_line", "chrom", "start", "strand", "strict_canonical"),
    col_types = readr::cols(biotype_folder = "c", cell_line = "c", chrom = "c",
                            start = "d", strand = "c", strict_canonical = "c")
  ) %>%
    dplyr::rename(biotype = biotype_folder) %>%
    dplyr::mutate(strict_canonical = toupper(trimws(strict_canonical)) == "TRUE")
  dup <- SITE_TABLE %>% dplyr::count(dplyr::across(dplyr::all_of(SITE_KEYS))) %>%
    dplyr::filter(n > 1)
  if (nrow(dup) > 0) stop("--site-table has ", nrow(dup), " duplicated site keys")
  message("  Site table: ", nrow(SITE_TABLE), " sites, ",
          sum(SITE_TABLE$strict_canonical), " strict_canonical == True")
}
SITE_MATCH <- list()

# Keep only replicate rows whose (biotype, cell_line, chrom, start, strand) is
# strict_canonical == True in the site table; record how unique sites matched.
apply_site_table <- function(rows, set_label) {
  if (is.null(SITE_TABLE)) return(rows)
  status <- rows %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(SITE_KEYS))) %>%
    dplyr::left_join(SITE_TABLE, by = SITE_KEYS) %>%
    dplyr::mutate(status = dplyr::case_when(
      is.na(strict_canonical) ~ "not_in_table",
      strict_canonical        ~ "strict_true",
      TRUE                    ~ "strict_false"))
  SITE_MATCH[[set_label]] <<- status %>%
    dplyr::count(biotype, cell_line, status) %>%
    tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    dplyr::mutate(set = set_label, .before = 1)
  keep <- dplyr::filter(status, status == "strict_true") %>%
    dplyr::select(dplyr::all_of(SITE_KEYS))
  out <- dplyr::semi_join(rows, keep, by = SITE_KEYS)
  message("  --site-table (", set_label, "): ", nrow(rows), " -> ", nrow(out),
          " replicate rows; unique sites: ",
          paste(names(table(status$status)), table(status$status), sep = " = ", collapse = ", "))
  out
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load sites (replicate rows collapsed to unique genomic sites)
# ══════════════════════════════════════════════════════════════════════════════

# Read every call file of `condition` for each biotype, apply `row_filter`,
# and return one row per (sample, site). Lists are bound explicitly once
# (fixes the original's `%>% bind_rows(Filter(Negate(is.null), .))`, which
# included every row twice).
load_call_rows <- function(calls_dir, biotypes, condition, row_filter) {
  per_biotype <- lapply(biotypes, function(bt) {
    bt_dir <- file.path(calls_dir, bt)
    if (!dir.exists(bt_dir)) return(NULL)

    # Classify every file by its real condition via get_condition() (the same
    # samples.tsv-backed lookup the rest of the pipeline uses).
    files <- list.files(bt_dir, pattern = "_taps\\.tsv$", full.names = TRUE)
    files <- files[vapply(files, function(f) {
      sample <- sub(paste0("_", bt, "_taps\\.tsv$"), "", basename(f))
      get_condition(sample) == condition
    }, logical(1))]
    if (length(files) == 0) return(NULL)

    per_file <- lapply(files, function(f) {
      sample <- sub(paste0("_", bt, "_taps\\.tsv$"), "", basename(f))
      tryCatch({
        df <- read_tsv(f, show_col_types = FALSE,
                       col_types = cols(chrom = col_character())) %>%
          row_filter() %>%
          dplyr::mutate(
            strand    = ifelse(!is.na(strand) & strand == "-", "-", "+"),
            biotype   = bt,
            sample    = sample,
            cell_line = get_cell_line(sample)
          ) %>%
          dplyr::select(dplyr::all_of(SITE_KEYS), sample)
        # Memory guard kept from the original (cap per file before
        # collapsing). No current call file has more than 20000 rows.
        if (nrow(df) > 20000) df <- dplyr::slice_sample(df, n = 20000)
        df
      }, error = function(e) {
        message("  WARNING: failed to read ", f, ": ", conditionMessage(e))
        NULL
      })
    })
    dplyr::bind_rows(Filter(Negate(is.null), per_file))
  })
  dplyr::bind_rows(Filter(Negate(is.null), per_biotype))
}

# Collapse replicate rows to unique sites. In pooled mode, cell_line is
# replaced by a single label first, so a site present in both lines counts once.
collapse_sites <- function(rows, pool) {
  if (pool) rows <- dplyr::mutate(rows, cell_line = POOLED_LABEL)
  rows %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(SITE_KEYS))) %>%
    dplyr::summarise(n_replicates = dplyr::n_distinct(sample), .groups = "drop")
}

message("Loading high-confidence (foreground) sites...")
message("  Filters: mod_rate >= ", opt$`min-mod`, " | padj < ", opt$`max-padj`,
        " | coverage >= ", opt$`min-cov`, " | snp_flag == PASS | condition: ",
        opt$condition)
fg_rows <- load_call_rows(
  CALLS_DIR, BIOTYPES_USE, opt$condition,
  function(df) dplyr::filter(df, snp_flag == "PASS", coverage >= opt$`min-cov`,
                             mod_rate >= opt$`min-mod`, padj < opt$`max-padj`)
)
if (is.null(fg_rows) || nrow(fg_rows) == 0) {
  stop("No high-confidence sites found. Check filters or that calling has completed.")
}
fg_rows <- apply_site_table(fg_rows, "foreground")
if (nrow(fg_rows) == 0) stop("No foreground sites left after --site-table filter.")

# `sites`: unique per cell line (used by the 06a/06b logos and, unless
# pooling, by 06c). `fg_sites`: the unit of the 06c enrichment panels.
sites    <- collapse_sites(fg_rows, pool = FALSE)
fg_sites <- collapse_sites(fg_rows, pool = POOL)

message("  Foreground replicate rows: ", nrow(fg_rows),
        " -> unique sites per cell line: ", nrow(sites),
        if (POOL) paste0(" -> pooled unique sites: ", nrow(fg_sites)) else "")
message("\nUnique foreground sites per biotype:")
print(as.data.frame(dplyr::count(fg_sites, biotype, cell_line)))


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Extract genomic sequences around each site
# ══════════════════════════════════════════════════════════════════════════════

# Map Ensembl-style chromosome names onto the BSgenome's UCSC seqnames.
#   MT / M                  -> chrM
#   already a seqname       -> unchanged (e.g. chr1, chrUn_KI270733v1)
#   1..22, X, Y             -> chr1 .. chrY
#   unplaced / unlocalized contigs, Ensembl "<ID>.<n>" (e.g. KI270733.1,
#   GL000220.1)             -> the unique seqname ending in "_<ID>v<n>",
#                              optionally followed by _random/_alt/_fix
#                              (e.g. chrUn_KI270733v1, chr1_KI270706v1_random),
#                              looked up in seqnames(genome) rather than
#                              assuming chrUn.
# Names that resolve to no seqname (or to more than one) return NA; callers
# must treat NA as an error, never drop such sites silently.
normalise_chrom <- function(chrom, genome = GENOME) {
  sn <- seqnames(genome)
  map_one <- function(x) {
    if (is.na(x)) return(NA_character_)
    if (x %in% c("MT", "M", "chrMT")) return("chrM")
    if (x %in% sn) return(x)
    if (paste0("chr", x) %in% sn) return(paste0("chr", x))
    if (grepl("^[A-Za-z]+[0-9]+\\.[0-9]+$", x)) {
      id  <- sub("\\.([0-9]+)$", "v\\1", x)
      hit <- grep(paste0("^chr[^_]+_", id, "(_random|_alt|_fix)?$"), sn, value = TRUE)
      if (length(hit) == 1) return(hit)
    }
    NA_character_
  }
  u <- unique(chrom)
  unname(setNames(vapply(u, map_one, character(1)), u)[chrom])
}

extract_sequences <- function(sites_df, genome, half_window, max_sites = opt$`max-sites`) {
  # half_window = 5 means we extract from pos-5 to pos+5 (11 nt total)
  # The modified C is at position half_window + 1 (1-based, centre)

  sites_df <- sites_df %>%
    dplyr::mutate(
      chrom_norm = normalise_chrom(chrom, genome),
      seq_start  = start - half_window,     # 0-based start
      seq_end    = start + half_window + 1  # 0-based end (exclusive)
    )

  # Never drop sites silently because their chromosome name is unknown.
  unmapped <- sites_df %>% dplyr::filter(is.na(chrom_norm))
  if (nrow(unmapped) > 0) {
    stop(nrow(unmapped), " site(s) on chromosome name(s) not found in the BSgenome: ",
         paste(utils::head(unique(unmapped$chrom), 20), collapse = ", "))
  }
  n_in       <- nrow(sites_df)
  sites_df   <- dplyr::filter(sites_df, seq_start >= 0)
  n_offstart <- n_in - nrow(sites_df)
  # Windows running past the contig end would make getSeq() fail for the
  # whole batch; drop and count them instead.
  n_pre_end  <- nrow(sites_df)
  sites_df   <- dplyr::filter(sites_df,
                              seq_end <= GenomeInfoDb::seqlengths(genome)[chrom_norm])
  n_offend   <- n_pre_end - nrow(sites_df)

  if (nrow(sites_df) == 0) return(character(0))

  # Limit to max_sites per call to avoid memory issues
  if (nrow(sites_df) > max_sites) {
    set.seed(42)
    sites_df <- dplyr::slice_sample(sites_df, n = max_sites)
    message("  Subsampled to ", max_sites, " sites")
  }

  message("  Extracting sequences for ", nrow(sites_df),
          " sites (window ±", half_window, " nt)...")

  strand_vec <- ifelse(!is.na(sites_df$strand) & sites_df$strand == "-", "-", "+")

  gr <- GRanges(
    seqnames = sites_df$chrom_norm,
    ranges   = IRanges(start = sites_df$seq_start + 1,
                       end   = sites_df$seq_start + (2 * half_window + 1)),
    strand   = strand_vec
  )

  # A getSeq() failure is an error, not an empty (silently skipped) panel.
  seqs <- as.character(getSeq(genome, gr))

  # Keep only sequences of correct length (some near chromosome ends may be shorter)
  expected_len <- 2 * half_window + 1
  n_fetched <- length(seqs)
  seqs      <- seqs[nchar(seqs) == expected_len]
  n_badlen  <- n_fetched - length(seqs)

  # Remove sequences with N
  n_beforeN <- length(seqs)
  seqs <- seqs[!grepl("N", seqs, ignore.case = TRUE)]
  n_withN   <- n_beforeN - length(seqs)
  seqs <- toupper(seqs)
  if (n_offstart + n_offend + n_badlen + n_withN > 0) {
    message("  Dropped: ", n_offstart, " window before contig start, ",
            n_offend, " past contig end, ", n_badlen, " truncated, ",
            n_withN, " containing N")
  }

  # BSgenome returns DNA (T); this is RNA — convert before anything downstream
  # touches these sequences, so every consumer sees the RNA alphabet.
  seqs <- gsub("T", "U", seqs, fixed = TRUE)

  message("  Valid sequences retained: ", length(seqs))
  seqs
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Build sequence logo per window per biotype
# ══════════════════════════════════════════════════════════════════════════════

# Manuscript-wide "Sequence" base palette (matches the RNA modification /
# biotype / sample / software legend set, not the AURORA report palette).
SEQ_BASE_COLOURS <- c(A = "#D7E199", C = "#F4A582", G = "#D6604D", U = "#8B5A7A")

dna_colours <- make_col_scheme(
  chars = names(SEQ_BASE_COLOURS),
  cols  = unname(SEQ_BASE_COLOURS)
)

make_logo_plot <- function(seqs, title, half_window, n_sites) {
  if (length(seqs) < 10) {
    message("  Too few sequences (", length(seqs), ") for logo — skipping")
    return(NULL)
  }

  total_width <- 2 * half_window + 1
  centre_pos  <- half_window + 1

  # x-axis labels: -hw ... -1, C, +1 ... +hw
  x_labels <- c(
    paste0("-", half_window:1),
    "C",
    paste0("+", 1:half_window)
  )

  p <- ggplot() +
    geom_logo(seqs,
              method    = "probability",
              seq_type  = "rna",
              col_scheme = dna_colours,
              font      = "roboto_medium") +
    annotate("rect",
             xmin = centre_pos - 0.5,
             xmax = centre_pos + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = SEQ_BASE_COLOURS[["C"]], alpha = 0.12) +
    annotate("segment",
             x = centre_pos, xend = centre_pos,
             y = 0, yend = 1,
             colour = SEQ_BASE_COLOURS[["C"]], linewidth = 0.4, linetype = "dashed") +
    scale_x_continuous(
      breaks = 1:total_width,
      labels = x_labels,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title    = title,
      subtitle = paste0("n = ", scales::comma(n_sites),
                        " unique high-confidence sites (mod_rate ≥ ",
                        opt$`min-mod`, ", padj < ", opt$`max-padj`, ")"),
      x        = paste0("Position relative to m5C (±", half_window, " nt)"),
      y        = "Base frequency"
    ) +
    theme_arial8 +
    theme(
      panel.grid.major.x = element_blank(),
      text          = element_text(size = 7, family = "Arial", colour = "black"),
      axis.text     = element_text(size = 7, family = "Arial", colour = "black"),
      axis.text.x   = element_text(size = 5.5, family = "Arial", colour = "black"),
      axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
      plot.title    = element_text(size = 7, family = "Arial", colour = "black"),
      plot.subtitle = element_text(size = 6, family = "Arial", colour = "grey50")
    )

  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Probability logos (06a, 06b) — unique sites, fixed ±10 nt
# ══════════════════════════════════════════════════════════════════════════════

for (hw in WINDOWS) {
  message("\n", strrep("=", 55))
  message("Window ±", hw, " nt")
  message(strrep("=", 55))

  # ── Combined logo (all biotypes) ───────────────────────────────────────────
  message("Combined (all biotypes):")
  seqs_all <- extract_sequences(sites, GENOME, hw)

  p_combined <- make_logo_plot(
    seqs  = seqs_all,
    title = paste0("Sequence context of m5C sites — all biotypes (±", hw, " nt)"),
    half_window = hw,
    n_sites     = length(seqs_all)
  )

  if (!is.null(p_combined)) {
    save_figure(
      p_combined,
      file.path(opt$figdir, "06a_sequence_context_all.pdf"),
      width = 5.2,
      height = 2.0
    )
    message("Figure 6a: combined sequence context — done")
  }

  # ── Per-biotype logos ─────────────────────────────────────────────────────
  bio_plots <- list()

  for (bt in BIOTYPES_USE) {
    bt_sites <- sites %>% dplyr::filter(biotype == bt)
    if (nrow(bt_sites) < 10) {
      message("  ", bt, ": too few sites (", nrow(bt_sites), ") — skipping")
      next
    }

    message(bt, " (", nrow(bt_sites), " sites):")
    seqs_bt <- extract_sequences(bt_sites, GENOME, hw)

    p_bt <- make_logo_plot(
      seqs        = seqs_bt,
      title       = paste0(bt, " m5C context (±", hw, " nt)"),
      half_window = hw,
      n_sites     = length(seqs_bt)
    )

    if (!is.null(p_bt)) bio_plots[[bt]] <- p_bt
  }

  if (length(bio_plots) >= 2) {
    panel <- wrap_plots(bio_plots, ncol = 1) +
      plot_annotation(
        title    = paste0("m5C sequence context per RNA biotype (±", hw, " nt)"),
        subtitle = paste0("Unique high-confidence sites: mod_rate ≥ ", opt$`min-mod`,
                          " | padj < ", opt$`max-padj`,
                          " | TET+PB condition"),
        theme    = theme_arial8
      )

    save_figure(
      panel,
      file.path(opt$figdir, "06b_sequence_context_by_biotype.pdf"),
      width  = 6.5,
      height = 2.0 * length(bio_plots)
    )
    message("Figure 6b: per-biotype sequence context — done")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Background-corrected enrichment (06c)
# ══════════════════════════════════════════════════════════════════════════════
# A flat-looking probability logo is only evidence of "no context bias" if it's
# compared against the base composition the assay would see anyway.
# Background = unique C sites from the SAME call files (same biotype,
# condition, coverage, snp filter) that are clearly NOT modified
# (mod_rate <= bg-max-mod), excluding any site that is foreground for the same
# panel (cell_line x biotype, or biotype when pooled).

message("\nBuilding background (unmodified C) unique-site set...")
bg_rows <- load_call_rows(
  CALLS_DIR, BIOTYPES_USE, opt$condition,
  function(df) dplyr::filter(df, snp_flag == "PASS", coverage >= opt$`min-cov`,
                             mod_rate <= opt$`bg-max-mod`)
)
bg_rows <- apply_site_table(bg_rows, "background")
bg_universe <- collapse_sites(bg_rows, pool = POOL) %>%
  dplyr::anti_join(fg_sites, by = SITE_KEYS)

# 3x oversample per panel (min 500, as in the original), from unique sites.
fg_panel_n <- dplyr::count(fg_sites, biotype, cell_line, name = "n_fg_sites")
set.seed(42)
bg_sites <- bg_universe %>%
  dplyr::inner_join(fg_panel_n, by = c("biotype", "cell_line")) %>%
  dplyr::group_by(biotype, cell_line) %>%
  dplyr::group_modify(~ {
    target <- max(3 * .x$n_fg_sites[1], 500)
    dplyr::slice_sample(.x, n = min(nrow(.x), target))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n_fg_sites)

message("  Background replicate rows: ", nrow(bg_rows),
        " -> unique, non-foreground sites: ", nrow(bg_universe),
        " -> sampled: ", nrow(bg_sites), " (vs ", nrow(fg_sites), " foreground)")

MC_B        <- 100000  # Monte Carlo tables per chi-squared test (min p = 1/(B+1))
hw_c        <- opt$`hw-c`
total_width <- 2 * hw_c + 1
centre_pos  <- hw_c + 1
offsets     <- seq(-hw_c, hw_c)
bases       <- c("A", "C", "G", "U")
# Label every other tick above ±5 nt so offset labels never collide.
x_labels    <- ifelse(offsets == 0, "C",
               ifelse(hw_c > 5 & offsets %% 2 != 0, "",
               ifelse(offsets > 0, paste0("+", offsets), as.character(offsets))))
PANEL_W_MM  <- if (hw_c <= 5) 70 else 85
PANEL_H_MM  <- 18
BASE_LETTER_COLOURS <- SEQ_BASE_COLOURS

position_base_counts <- function(seqs) {
  m <- do.call(rbind, strsplit(seqs, ""))
  per_position <- lapply(seq_len(ncol(m)), function(j) {
    tab <- table(factor(m[, j], levels = bases))
    data.frame(position = j, base = bases, count = as.integer(tab))
  })
  dplyr::bind_rows(per_position)
}

compute_enrich <- function(seqs_fg, seqs_bg) {
  fg_bc <- position_base_counts(seqs_fg) %>%
    dplyr::mutate(freq = count / length(seqs_fg) * 100)
  bg_bc <- position_base_counts(seqs_bg) %>%
    dplyr::mutate(freq = count / length(seqs_bg) * 100)

  # Per-position chi-squared on the 2x4 (fg/bg x A/C/G/U) table with a
  # Monte Carlo p-value (B = MC_B tables, conditional on the margins), so the
  # test stays valid when expected counts are small; BH-adjusted across the
  # window. The seed is set once per panel, before the position loop, so
  # results are reproducible. A base absent from BOTH sets at a position (an
  # all-zero column) carries no information and is dropped before testing;
  # positions with fewer than two observed bases (the all-C centre) are NA and
  # are not counted as tests.
  set.seed(42)
  pos_pvals <- lapply(seq_len(total_width), function(j) {
    fg_j <- fg_bc$count[fg_bc$position == j]
    bg_j <- bg_bc$count[bg_bc$position == j]
    tab  <- rbind(fg_j, bg_j)
    tab  <- tab[, colSums(tab) > 0, drop = FALSE]
    ct <- if (ncol(tab) < 2) NULL else
      stats::chisq.test(tab, simulate.p.value = TRUE, B = MC_B)
    data.frame(position     = j,
               statistic    = if (is.null(ct)) NA_real_ else unname(ct$statistic),
               pval         = if (is.null(ct)) NA_real_ else ct$p.value,
               min_expected = if (is.null(ct)) NA_real_ else min(ct$expected))
  }) %>% dplyr::bind_rows() %>%
    dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"),
                  sig  = !is.na(padj) & padj < 0.05)

  enrich_all <- fg_bc %>%
    dplyr::select(position, base, fg_count = count, fg_freq = freq) %>%
    dplyr::left_join(dplyr::select(bg_bc, position, base, bg_count = count, bg_freq = freq),
                     by = c("position", "base")) %>%
    dplyr::mutate(enrichment = fg_freq - bg_freq)
  enrich <- dplyr::filter(enrich_all, position != centre_pos)  # skip tautological C

  sig_y <- enrich %>%
    dplyr::group_by(position) %>%
    dplyr::summarise(ytop = sum(pmax(enrichment, 0)) + 0.6, .groups = "drop") %>%
    dplyr::inner_join(dplyr::filter(pos_pvals, sig), by = "position")

  list(enrich = enrich, enrich_all = enrich_all, sig_y = sig_y, pos_pvals = pos_pvals)
}

# True two-sample logo (Vacic et al. 2006 style): real ggseqlogo letter-glyph
# polygons on both sides of zero; enriched and depleted bases are rendered as
# two separate probability-normalised logos rescaled to percentage points.
render_diff_plot <- function(d, title, n_fg, n_bg) {
  real_positions <- setdiff(1:total_width, centre_pos)

  wide <- d$enrich %>%
    dplyr::mutate(pos_val = pmax(enrichment, 0), neg_val = pmax(-enrichment, 0))

  mk_matrix <- function(value_col) {
    w <- tidyr::pivot_wider(wide, id_cols = base, names_from = position,
                            values_from = dplyr::all_of(value_col), values_fill = 0)
    m <- as.matrix(w[, as.character(real_positions)])
    rownames(m) <- w$base
    m[bases, , drop = FALSE]
  }
  mat_pos <- mk_matrix("pos_val")
  mat_neg <- mk_matrix("neg_val")

  total_pos <- setNames(colSums(mat_pos), real_positions)
  total_neg <- setNames(colSums(mat_neg), real_positions)

  # Depleted letters stay upright and sit below zero ((y-1)*mag).
  build_half <- function(mat, totals, upright_below, half_tag) {
    poly <- geom_logo(mat, method = "probability", seq_type = "rna", plot = FALSE)
    col_to_realpos <- setNames(real_positions, seq_along(real_positions))
    poly %>%
      dplyr::mutate(
        real_pos = col_to_realpos[as.character(position)],
        x        = x - position + real_pos,
        y        = if (upright_below) (y - 1) * totals[as.character(real_pos)]
                   else y * totals[as.character(real_pos)],
        half     = half_tag
      )
  }
  ref_label <- data.frame(x = centre_pos, y = 0, label = "m5C")

  polys <- dplyr::bind_rows(
    build_half(mat_pos, total_pos, FALSE, "pos"),
    build_half(mat_neg, total_neg, TRUE, "neg")
  ) %>%
    # glyphs with holes are one self-intersecting path; row order is load-bearing
    dplyr::arrange(real_pos, letter, half, order)

  ggplot(polys, aes(x = x, y = y, group = interaction(real_pos, letter, half),
                    fill = letter)) +
    geom_polygon(rule = "evenodd") +
    geom_text(data = ref_label, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2.46, fontface = "bold",
              colour = SEQ_BASE_COLOURS[["C"]], family = "Arial") +
    geom_text(data = d$sig_y, aes(x = position, y = ytop, label = "*"),
              inherit.aes = FALSE, size = 2.3, colour = "black", family = "Arial") +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = centre_pos, colour = "grey80", linewidth = 0.3,
              linetype = "dashed") +
    scale_fill_manual(values = BASE_LETTER_COLOURS, guide = "none") +
    scale_x_continuous(breaks = 1:total_width, labels = x_labels, expand = c(0.01, 0.01)) +
    labs(
      title = title,
      subtitle = paste0("n = ", scales::comma(n_fg), " m5C vs ",
                        scales::comma(n_bg), " unmodified C (unique sites)",
                        " | * BH Monte Carlo χ² padj < 0.05"),
      # short y title: the rotated label must fit within the panel height,
      # otherwise it overruns the subtitle above it
      x = "Distance from m5C (nt)", y = "Δ freq. (pp)"
    ) +
    theme_arial8 +
    theme(
      panel.grid.major.x    = element_blank(),
      plot.title.position   = "plot",
      text          = element_text(size = 7, family = "Arial", colour = "black"),
      axis.text     = element_text(size = 7, family = "Arial", colour = "black"),
      axis.title    = element_text(size = 7, family = "Arial", colour = "black", face = "bold"),
      plot.title    = element_text(size = 7, family = "Arial", colour = "black"),
      plot.subtitle = element_text(size = 6, family = "Arial", colour = "grey50")
    )
}

PANEL_GROUPS <- if (POOL) POOLED_LABEL else c("Caco2", "HEK")

diff_data <- list()
for (cl in PANEL_GROUPS) {
  for (bt in BIOTYPES_USE) {
    bt_fg_sites <- fg_sites %>% dplyr::filter(biotype == bt, cell_line == cl)
    bt_bg_sites <- bg_sites %>% dplyr::filter(biotype == bt, cell_line == cl)
    key <- paste(bt, cl, sep = " | ")
    if (nrow(bt_fg_sites) < 10 || nrow(bt_bg_sites) < 10) {
      message("  ", key, ": too few fg/bg sites for enrichment — skipping")
      next
    }
    message(key, ": ", nrow(bt_fg_sites), " fg / ", nrow(bt_bg_sites), " bg unique sites")
    seqs_bt_fg <- extract_sequences(bt_fg_sites, GENOME, hw_c, max_sites = nrow(bt_fg_sites))
    seqs_bt_bg <- extract_sequences(bt_bg_sites, GENOME, hw_c, max_sites = nrow(bt_bg_sites))
    if (length(seqs_bt_fg) < 10 || length(seqs_bt_bg) < 10) next

    diff_data[[key]] <- c(
      compute_enrich(seqs_bt_fg, seqs_bt_bg),
      list(n_fg = length(seqs_bt_fg), n_bg = length(seqs_bt_bg),
           n_fg_sites = nrow(bt_fg_sites), n_bg_sites = nrow(bt_bg_sites),
           biotype = bt, cell_line = cl)
    )
  }
}

if (length(diff_data) == 0) stop("No 06c panels could be computed.")

# ── Tables (always written) ──────────────────────────────────────────────────
sig_table <- lapply(diff_data, function(d) {
  d$pos_pvals %>%
    dplyr::mutate(biotype = d$biotype, cell_line = d$cell_line, hw_c = hw_c,
                  offset = position - centre_pos,
                  test = paste0("chisq_monte_carlo_B", MC_B, "_seed42"),
                  n_fg = d$n_fg, n_bg = d$n_bg)
}) %>% dplyr::bind_rows() %>%
  dplyr::select(biotype, cell_line, hw_c, offset, position, test, statistic, pval, padj,
                sig, min_expected, n_fg, n_bg) %>%
  dplyr::arrange(cell_line, biotype, position)
sig_out <- file.path(opt$figdir, "06c_sequence_context_enrichment_significance.tsv")
readr::write_tsv(sig_table, sig_out)
message("  Saved: ", sig_out)

freq_table <- lapply(diff_data, function(d) {
  d$enrich_all %>%
    dplyr::mutate(biotype = d$biotype, cell_line = d$cell_line, hw_c = hw_c,
                  offset = position - centre_pos)
}) %>% dplyr::bind_rows() %>%
  dplyr::select(biotype, cell_line, hw_c, offset, position, base,
                fg_count, bg_count, fg_freq, bg_freq, enrichment) %>%
  dplyr::arrange(cell_line, biotype, position, base)
freq_out <- file.path(opt$figdir, "06c_sequence_context_enrichment_frequencies.tsv")
readr::write_tsv(freq_table, freq_out)
message("  Saved: ", freq_out)

# ── Figure ───────────────────────────────────────────────────────────────────
fixed_gtable <- function(p) {
  gt <- ggplotGrob(p)
  is_panel  <- grepl("^panel($|-)", gt$layout$name)
  panel_row <- unique(gt$layout$t[is_panel])
  panel_col <- unique(gt$layout$l[is_panel])
  gt$widths[panel_col]  <- unit(PANEL_W_MM, "mm")
  gt$heights[panel_row] <- unit(PANEL_H_MM, "mm")
  gt
}
gt_size_in <- function(gt) {
  pdf(NULL)
  on.exit(dev.off())
  c(w = sum(grid::convertWidth(gt$widths, "in", valueOnly = TRUE)),
    h = sum(grid::convertHeight(gt$heights, "in", valueOnly = TRUE)))
}

panel_title <- function(d) {
  if (POOL) paste0(d$biotype, " · ", POOLED_LABEL, " pooled (±", hw_c, " nt)")
  else      paste0(d$biotype, " · ", d$cell_line, " (±", hw_c, " nt)")
}
diff_gts <- lapply(diff_data, function(d) {
  fixed_gtable(render_diff_plot(d, panel_title(d), d$n_fg, d$n_bg))
})

stack_column <- function(gts) do.call(rbind, c(unname(gts), list(size = "max")))

two_column <- !POOL && hw_c > 5
if (two_column) {
  cols <- lapply(PANEL_GROUPS, function(cl) {
    keys <- paste(BIOTYPES_USE, cl, sep = " | ")
    keys <- keys[keys %in% names(diff_gts)]
    if (length(keys) == 0) NULL else stack_column(diff_gts[keys])
  })
  cols  <- Filter(Negate(is.null), cols)
  sizes <- vapply(cols, gt_size_in, numeric(2))
  total_w_in <- sum(sizes["w", ])
  total_h_in <- max(sizes["h", ])
  # Pad shorter columns at the bottom so panels stay top-aligned.
  figure <- gridExtra::arrangeGrob(
    grobs  = lapply(seq_along(cols), function(i) {
      gridExtra::arrangeGrob(cols[[i]], grid::nullGrob(),
                             heights = unit(c(sizes["h", i], total_h_in - sizes["h", i]), "in"))
    }),
    ncol   = length(cols),
    widths = unit(sizes["w", ], "in")
  )
} else {
  figure <- stack_column(diff_gts)
  sz <- gt_size_in(figure)
  total_w_in <- sz[["w"]]
  total_h_in <- sz[["h"]]
}

base_out <- file.path(opt$figdir, "06c_sequence_context_enrichment")
ggsave(paste0(base_out, ".pdf"), figure, width = total_w_in, height = total_h_in,
       device = cairo_pdf, family = "Arial")
ggsave(paste0(base_out, ".png"), figure, width = total_w_in, height = total_h_in, dpi = 300)
ggsave(paste0(base_out, ".svg"), figure, width = total_w_in, height = total_h_in, device = "svg")
message("  Saved: ", base_out, ".{pdf,png,svg} | ", length(diff_gts), " panels, ",
        if (two_column) "two-column" else "one column", ", panel ",
        PANEL_W_MM, " x ", PANEL_H_MM, " mm, ",
        round(total_w_in * 25.4), " x ", round(total_h_in * 25.4), " mm overall")

# ── Summary table of unique sites used ───────────────────────────────────────
summary_out <- file.path(opt$figdir, "06_sequence_context_site_summary.tsv")
dplyr::bind_rows(lapply(diff_data, function(d) data.frame(
  biotype = d$biotype, cell_line = d$cell_line, hw_c = hw_c,
  fg_unique_sites = d$n_fg_sites, bg_unique_sites_sampled = d$n_bg_sites,
  fg_sequences = d$n_fg, bg_sequences = d$n_bg
))) %>%
  dplyr::left_join(
    dplyr::count(bg_universe, biotype, cell_line, name = "bg_unique_sites_available"),
    by = c("biotype", "cell_line")
  ) %>%
  readr::write_tsv(summary_out)
message("\nSite summary: ", summary_out)

if (length(SITE_MATCH) > 0) {
  match_out <- file.path(opt$figdir, "06_sequence_context_site_table_match.tsv")
  dplyr::bind_rows(SITE_MATCH) %>%
    dplyr::arrange(set, biotype, cell_line) %>%
    readr::write_tsv(match_out)
  message("Site-table match summary (unique sites per cell line, before pooling): ", match_out)
}
message("Figure 6 sequence-context generation complete.")
