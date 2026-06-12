# -*- coding: utf-8 -*-
"""
09_compare.py
Benchmarking comparison: custom TAPS pipeline vs rastair, asTair, Bismark.

Tool characteristics:
    rastair  — TAPS-native, CpG-only, Bowtie1 BAMs. Direct caller comparison.
               Comparison restricted to CpG-context sites from custom pipeline.

    asTair   — TAPS-native, all C contexts, Bowtie1 biotype BAMs. Fully matched
               aligner + input. Most direct like-for-like comparison.

    Bismark  — Bisulfite caller on Bowtie1 alignments. Included as contrast.
               CHEMISTRY INVERSION: Bismark reads C→T as UNMETHYLATED.
               TAPS reads C→T as METHYLATED. They are chemical inverses.
               Before any comparison, Bismark mod_rate is inverted:
                   bismark_taps_equiv = 1.0 - bismark_methylation_rate
               After inversion, concordance and correlation are computed
               as normal. Anti-correlation before inversion is expected and
               serves as a positive control that the chemistry is correct.

Metrics per tool per condition per biotype:
    1. Site-level concordance : Jaccard, recall, shared site count
    2. Mod-rate correlation   : Pearson r, Spearman rho at shared sites
    3. Pre/post inversion     : Bismark reported both ways for transparency

Output: 09.compare/
    concordance_summary.tsv
    correlation_summary.tsv
    shared_sites_<condition>_<biotype>_<tool>.tsv

Usage:
    python 09_compare.py
    python 09_compare.py --condition treat --biotype miRNA --tool astair
"""

import argparse
import glob
import logging
import os
import subprocess
from pathlib import Path

import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
CALLS_DIR = "07.taps_calls"
BENCH_DIR = "08.benchmark"
OUT_DIR   = "09.compare"

CONDITIONS = [
    "no-treat_Ctrl_Caco2", "no-treat_Ctrl_HEK",
    "pb_Ctrl_Caco2",        "pb_Ctrl_HEK",
    "treat_Caco2",          "treat_HEK",
]
BIOTYPES = ["rRNA", "miRNA", "tRNA", "snoRNA", "snRNA", "piRNA", "lncRNA", "other"]
TOOLS    = ["rastair", "rastair_all", "astair", "bismark"]

# Bismark chemistry is inverted vs TAPS — applied before comparison
INVERT_CHEMISTRY = {"bismark": True, "rastair": False, "rastair_all": False, "astair": False}

# rastair is CpG-only — comparison restricted to CpG context sites
CpG_ONLY = {"rastair": True, "rastair_all": False, "astair": False, "bismark": False}


# ════════════════════════════════════════════════════════════════════════════
# Loaders
# ════════════════════════════════════════════════════════════════════════════

def load_custom(condition, biotype):
    """Load all replicates for a condition+biotype, PASS sites only."""
    pattern = f"{CALLS_DIR}/{biotype}/{condition}*_{biotype}_taps.tsv"
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        dfs.append(df)
    if not dfs:
        return None
    combined = pd.concat(dfs, ignore_index=True)
    # Keep PASS sites only
    if "snp_flag" in combined.columns:
        combined = combined[combined["snp_flag"] == "PASS"].copy()
    combined["site_key"] = (
        combined["chrom"].astype(str) + ":" +
        combined["start"].astype(str)
    )
    return combined


def load_rastair(condition, subdir="cpg"):
    """
    Load rastair 2.1.1 BED output.
    subdir: "cpg" for CpG-only run, "all" for all-context run.
    Actual columns (with header line starting #):
        #chr, start, end, name, beta_est, strand, unmod, mod,
        no_snp, snp, coverage, genotype, gt_p_score, gt_conf_score, cpg
    mod_rate = beta_est (0.0-1.0)
    coverage = coverage column
    """
    pattern = f"{BENCH_DIR}/rastair/{subdir}/{condition}*/*.bed.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        pattern = f"{BENCH_DIR}/rastair/{subdir}/{condition}*/*.bed"
        files = sorted(glob.glob(pattern))
    # fallback to old flat structure
    if not files:
        pattern = f"{BENCH_DIR}/rastair/{condition}*/*.bed.gz"
        files = sorted(glob.glob(pattern))
    if not files:
        pattern = f"{BENCH_DIR}/rastair/{condition}*/*.bed"
        files = sorted(glob.glob(pattern))
    if not files:
        return None
    dfs = []
    for f in files:
        compression = "gzip" if f.endswith(".gz") else None
        try:
            df = pd.read_csv(f, sep="\t", compression=compression,
                             dtype={"#chr": str})
            # Strip # from column name
            df.columns = [c.lstrip("#") for c in df.columns]
            df = df.rename(columns={
                "chr":      "chrom",
                "beta_est": "mod_rate",
            })
            # mod_rate is already 0-1
            df["mod_rate"] = pd.to_numeric(df["mod_rate"], errors="coerce")
            df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
            df = df.dropna(subset=["mod_rate", "coverage"])
            dfs.append(df[["chrom", "start", "mod_rate", "coverage"]])
        except Exception as e:
            log.warning("rastair load error %s: %s", f, e)
            continue
    if not dfs:
        return None
    combined = pd.concat(dfs, ignore_index=True)
    combined["chrom"] = combined["chrom"].astype(str)
    combined["site_key"] = combined["chrom"] + ":" + combined["start"].astype(str)
    return combined[["site_key", "chrom", "start", "mod_rate", "coverage"]].copy()


def load_astair(condition, biotype, context="all"):
    """
    Load asTair mCaller output for a specific biotype and context.
    """
    # asTair names output from BAM basename: <bam_name>_<context>.mCaller.gz
    # BAM basename format: <condition>_<biotype>.sorted
    # Output: <condition>_<biotype>.sorted_<context>.mCaller.gz
    pattern = f"{BENCH_DIR}/astair/{biotype}/{condition}*_{biotype}.sorted_mCtoT_{context}.mods.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    dfs = []
    for f in files:
        # asTair column names: #CHROM START END MOD_LEVEL MOD UNMOD TOTAL_DEPTH
        df = pd.read_csv(f, sep="\t", compression="gzip",
                         comment=None, header=0)
        # Strip leading # from column names (asTair uses #CHROM)
        df.columns = [str(c).lstrip("#") for c in df.columns]
        # Rename to standard names
        rename = {
            "CHROM": "chrom", "START": "start", "END": "end",
            "MOD_LEVEL": "mod_rate", "MOD": "mod_count",
            "UNMOD": "unmod_count", "TOTAL_DEPTH": "coverage",
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
        if "mod_rate" not in df.columns:
            continue  # skip if rename failed
        # Drop zero-coverage sites (MOD_LEVEL == "*")
        df = df[df["mod_rate"].astype(str) != "*"].copy()
        df["mod_rate"] = pd.to_numeric(df["mod_rate"], errors="coerce")
        df["chrom"]    = df["chrom"].astype(str)
        df = df.dropna(subset=["mod_rate"])
        dfs.append(df)
    if not dfs:
        return None
    combined = pd.concat(dfs, ignore_index=True)
    combined["site_key"] = combined["chrom"].astype(str) + ":" + combined["start"].astype(str)
    return combined[["site_key", "chrom", "start", "mod_rate", "coverage"]].copy()


def load_bismark(condition, biotype):
    """
    Load Bismark methylation extractor output (CX_report or bismark.cov).
    Handles both --CX_context (all contexts) and CpG-only output.

    IMPORTANT: Bismark mod_rate is inverted before comparison.
    bismark_methylation = fraction of reads showing C (unmodified in bisulfite)
    In TAPS, C = unmodified, T = modified.
    Therefore: taps_equivalent = 1.0 - bismark_methylation_rate
    """
    # Try CX context report first (produced by --CX_context flag)
    pattern = f"{BENCH_DIR}/bismark/{condition}*/*.CX_report.txt.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        # Fall back to bedGraph coverage files
        pattern = f"{BENCH_DIR}/bismark/{condition}*/*.bismark.cov.gz"
        files = sorted(glob.glob(pattern))
    if not files:
        return None

    dfs = []
    for f in files:
        if f.endswith(".CX_report.txt.gz"):
            # CX report format: chrom, pos, strand, count_M, count_U, context, trinucleotide
            df = pd.read_csv(f, sep="\t", compression="gzip", header=None,
                             names=["chrom","start","strand","count_M","count_U","context","trinuc"])
            df["coverage"] = df["count_M"] + df["count_U"]
            df = df[df["coverage"] > 0].copy()
            df["bismark_rate"] = df["count_M"] / df["coverage"]
        else:
            # bismark.cov format: chrom, start, end, methylation_pct, count_M, count_U
            df = pd.read_csv(f, sep="\t", compression="gzip", header=None,
                             names=["chrom","start","end","meth_pct","count_M","count_U"])
            df["coverage"] = df["count_M"] + df["count_U"]
            df["bismark_rate"] = df["meth_pct"] / 100.0
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)

    # ── Chemistry inversion ───────────────────────────────────────────────────
    # Bismark bismark_rate = fraction METHYLATED in bisulfite logic
    #                      = fraction of C reads (unmodified C in bisulfite)
    # In TAPS:  C = unmodified,  T = modified
    # In BS-seq: C = methylated, T = unmethylated  (opposite!)
    # Therefore TAPS mod_rate ≈ 1.0 - bismark_rate at true modification sites
    combined["mod_rate_raw"]     = combined["bismark_rate"]          # before inversion
    combined["mod_rate"]         = 1.0 - combined["bismark_rate"]    # after inversion
    combined["chemistry_note"]   = "bismark_inverted_1_minus_rate"

    combined["site_key"] = combined["chrom"].astype(str) + ":" + combined["start"].astype(str)
    return combined[["site_key", "chrom", "start", "mod_rate", "mod_rate_raw", "coverage"]].copy()


# ════════════════════════════════════════════════════════════════════════════
# Metrics
# ════════════════════════════════════════════════════════════════════════════

def concordance(df_a, df_b, name_a, name_b, condition, biotype, tool):
    sa = set(df_a["site_key"].dropna())
    sb = set(df_b["site_key"].dropna())
    shared = len(sa & sb)
    union  = len(sa | sb)
    return {
        "condition":       condition,
        "biotype":         biotype,
        "tool":            tool,
        "sites_custom":    len(sa),
        "sites_tool":      len(sb),
        "shared":          shared,
        "union":           union,
        "jaccard":         round(shared / union, 4) if union else 0.0,
        "recall_custom":   round(shared / len(sa), 4) if sa else 0.0,
        "recall_tool":     round(shared / len(sb), 4) if sb else 0.0,
        "chemistry_inverted": INVERT_CHEMISTRY.get(tool, False),
    }


def correlation(df_a, df_b, name_a, name_b, condition, biotype, tool):
    merged = df_a[["site_key", "mod_rate"]].merge(
        df_b[["site_key", "mod_rate"]], on="site_key",
        suffixes=("_custom", f"_{name_b}")
    ).dropna()

    result = {
        "condition":          condition,
        "biotype":            biotype,
        "tool":               tool,
        "n_shared":           len(merged),
        "chemistry_inverted": INVERT_CHEMISTRY.get(tool, False),
    }

    if len(merged) >= 5:
        a = merged["mod_rate_custom"]
        b = merged[f"mod_rate_{name_b}"]
        pr, pp = stats.pearsonr(a, b)
        sr, sp = stats.spearmanr(a, b)
        result.update({
            "pearson_r":    round(pr, 4),
            "pearson_p":    round(pp, 6),
            "spearman_rho": round(sr, 4),
            "spearman_p":   round(sp, 6),
        })
    else:
        result.update({
            "pearson_r": None, "pearson_p": None,
            "spearman_rho": None, "spearman_p": None,
        })
    return result, merged


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--condition", default="all")
    p.add_argument("--biotype",   default="all")
    p.add_argument("--tool",      default="all")
    return p.parse_args()


def main():
    args = parse_args()
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    conditions = CONDITIONS if args.condition == "all" else [args.condition]
    biotypes   = BIOTYPES   if args.biotype   == "all" else [args.biotype]
    tools      = TOOLS      if args.tool      == "all" else [args.tool]

    conc_rows = []
    corr_rows = []

    for condition in conditions:
        for biotype in biotypes:

            custom = load_custom(condition, biotype)
            if custom is None:
                log.info("No custom calls: %s / %s — skipping", condition, biotype)
                continue

            log.info("%s / %s : %d PASS sites", condition, biotype, len(custom))

            for tool in tools:

                # ── Load tool output ──────────────────────────────────────────
                tool_df = None
                if tool == "rastair":
                    tool_df = load_rastair(condition, subdir="cpg")
                    # rastair is CpG-only: restrict custom to CpG sites
                    if tool_df is not None and "context" in custom.columns:
                        custom_sub = custom[custom["context"].str[1:3] == "CG"].copy()
                    else:
                        custom_sub = custom
                elif tool == "rastair_all":
                    tool_df = load_rastair(condition, subdir="all")
                    # all-context: compare against full custom call set
                    custom_sub = custom
                elif tool == "astair":
                    tool_df    = load_astair(condition, biotype, context="all")
                    custom_sub = custom
                elif tool == "bismark":
                    tool_df    = load_bismark(condition, biotype)
                    custom_sub = custom

                if tool_df is None:
                    log.info("  %-10s : output not found — skipping", tool)
                    continue

                log.info("  %-10s : %d sites%s",
                         tool, len(tool_df),
                         " [chemistry inverted]" if INVERT_CHEMISTRY[tool] else "")

                # ── Concordance ───────────────────────────────────────────────
                c = concordance(custom_sub, tool_df,
                                "custom", tool,
                                condition, biotype, tool)
                conc_rows.append(c)
                log.info("             shared=%d  Jaccard=%.4f",
                         c["shared"], c["jaccard"])

                # ── Correlation ───────────────────────────────────────────────
                r, merged = correlation(custom_sub, tool_df,
                                        "custom", tool,
                                        condition, biotype, tool)
                corr_rows.append(r)
                log.info("             Pearson=%-6s  Spearman=%-6s  n=%d",
                         r.get("pearson_r", "NA"),
                         r.get("spearman_rho", "NA"),
                         r["n_shared"])

                # Save shared sites table
                if len(merged) > 0:
                    out_f = f"{OUT_DIR}/shared_{condition}_{biotype}_{tool}.tsv"
                    merged.to_csv(out_f, sep="\t", index=False)

    # ── Write summaries ───────────────────────────────────────────────────────
    if conc_rows:
        out = f"{OUT_DIR}/concordance_summary.tsv"
        pd.DataFrame(conc_rows).to_csv(out, sep="\t", index=False)
        log.info("Concordance → %s", out)

    if corr_rows:
        out = f"{OUT_DIR}/correlation_summary.tsv"
        pd.DataFrame(corr_rows).to_csv(out, sep="\t", index=False)
        log.info("Correlation → %s", out)

    # ── MultiQC ───────────────────────────────────────────────────────────────
    mq_dirs = " ".join([BENCH_DIR, CALLS_DIR, "02.fastqc", "03.trimGalore"])
    mq_out  = f"{OUT_DIR}/multiqc_compare"
    os.makedirs(mq_out, exist_ok=True)
    cmd = (f"multiqc {mq_dirs} --outdir {mq_out} "
           f"--filename multiqc_pipeline_compare "
           f"--title 'akschneider TAPS-RNA: rastair vs asTair vs Bismark'")
    log.info("Running MultiQC...")
    subprocess.run(cmd, shell=True, check=False)

    log.info("Done → %s/", OUT_DIR)


if __name__ == "__main__":
    main()
