# -*- coding: utf-8 -*-
"""
07_taps_calling.py
TAPS m5C caller — chromosome-parallel with SNP filtering BEFORE calling.

SNP filtering order (correct):
    At each reference-C position, check SNP flags FIRST.
    If flagged → skip position entirely (no counts, no test, no output).
    Binomial test + BH correction run only on clean PASS positions.

    This matters because:
      a) SNP sites have mod_rate ~0.50 (het) or ~1.0 (hom) and produce
         highly significant p-values that pollute the BH correction pool,
         inflating or deflating adjusted p-values for real m5C sites.
      b) Running the test on SNPs wastes compute and produces misleading output.

Three-layer SNP filtering:

    Layer 1 — dbSNP (common population variants, AF >= 0.01):
        Known C→T or G→A SNPs from dbSNP hg38.
        Flagged: SNP_KNOWN

    Layer 2 — Cell-line-specific SNPs (from no-treat BAMs):
        C→T or G→A sites at AF >= 0.20 in the untreated condition.
        HEK293 and Caco2 carry mutations not in population dbSNP.
        Flagged: SNP_SAMPLE

    Layer 3 — Heterozygosity flag (from no-treat BAM rates):
        Sites where no-treat mod_rate >= het_threshold (default 0.40).
        Catches het SNPs not captured by AF threshold alone.
        Flagged: SNP_HET

    Sites with multiple flags: SNP_MULTI
    Sites with no flag: PASS → these are the only sites that proceed to calling.

Output columns:
    chrom, start, end, context, mod_count, unmod_count, coverage,
    mod_rate, pvalue, padj, snp_flag

Usage:
    python 07_taps_calling.py \
        --bam          sample_rRNA.sorted.bam \
        --fasta        hg38.fa \
        --out          sample_rRNA_taps.tsv \
        --dbsnp-vcf    06.snp_resources/dbsnp_hg38_CT_GA_af01.vcf.gz \
        --sample-snp-bed 06.snp_resources/sample_snps_HEK.bed \
        --cell-line    HEK \
        --min-cov      10 \
        --threads      4
"""

import argparse
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pysam
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ── SNP flag constants ────────────────────────────────────────────────────────
PASS       = "PASS"
SNP_KNOWN  = "SNP_KNOWN"
SNP_SAMPLE = "SNP_SAMPLE"
SNP_HET    = "SNP_HET"
SNP_MULTI  = "SNP_MULTI"


# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam",             required=True)
    p.add_argument("--fasta",           required=True)
    p.add_argument("--out",             required=True)
    p.add_argument("--min-qual",        type=int,   default=20)
    p.add_argument("--min-mapq",        type=int,   default=10)
    p.add_argument("--min-cov",         type=int,   default=5)
    p.add_argument("--context",         default="ALL",
                   choices=["ALL", "CpG", "CHH", "CHG"])
    p.add_argument("--threads",         type=int,   default=4)
    p.add_argument("--dbsnp-vcf",       default=None)
    p.add_argument("--sample-snp-bed",  default=None)
    p.add_argument("--het-threshold",   type=float, default=0.40)
    p.add_argument("--cell-line",       default="unknown")
    p.add_argument("--background-rate", type=float, default=0.005,
                   help="Expected C→T rate for unmodified C in binomial null (default: 0.005)")
    return p.parse_args()


# ════════════════════════════════════════════════════════════════════════════
# SNP resource loaders
# ════════════════════════════════════════════════════════════════════════════

def load_dbsnp(vcf_path):
    """
    Load known C→T / G→A positions from tabix-indexed dbSNP VCF.
    Returns frozenset of 'chrom:pos0' strings (0-based).
    frozenset is hashable, immutable, and pickles efficiently for workers.
    """
    if not vcf_path or not Path(vcf_path).exists():
        print(f"[snp] dbSNP VCF not found: {vcf_path} — Layer 1 disabled", flush=True)
        return frozenset()
    positions = set()
    try:
        vcf = pysam.VariantFile(vcf_path)
        for rec in vcf.fetch():
            ref  = rec.ref
            alts = [str(a) for a in rec.alts] if rec.alts else []
            if (ref == "C" and "T" in alts) or (ref == "G" and "A" in alts):
                positions.add(f"{rec.chrom}:{rec.pos - 1}")  # VCF 1-based → 0-based
        vcf.close()
    except Exception as e:
        print(f"[snp] WARNING: Could not read dbSNP VCF: {e}", flush=True)
        return frozenset()
    print(f"[snp] Layer 1 (dbSNP):      {len(positions):,} positions loaded", flush=True)
    return frozenset(positions)


def load_sample_snp_bed(bed_path, het_threshold):
    """
    Load cell-line-specific SNP BED.
    Returns two frozensets:
        sample_snps   — all positions in the BED (Layer 2, AF >= min_af at build time)
        het_positions — positions where AF >= het_threshold (Layer 3)
    """
    if not bed_path or not Path(bed_path).exists():
        print(f"[snp] Sample SNP BED not found: {bed_path} — Layers 2+3 disabled", flush=True)
        return frozenset(), frozenset()

    sample_snps   = set()
    het_positions = set()

    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            key = f"{parts[0]}:{parts[1]}"
            sample_snps.add(key)
            # Column 5 (0-indexed) is allele_freq written by 06_build_snp_blacklist.py
            if len(parts) >= 6:
                try:
                    if float(parts[5]) >= het_threshold:
                        het_positions.add(key)
                except ValueError:
                    pass

    print(f"[snp] Layer 2 (sample SNPs): {len(sample_snps):,} positions loaded", flush=True)
    print(f"[snp] Layer 3 (het≥{het_threshold}):   {len(het_positions):,} positions loaded", flush=True)
    return frozenset(sample_snps), frozenset(het_positions)


def assign_snp_flag(key, dbsnp, sample_snps, het_positions):
    """
    Assign SNP flag for a position key 'chrom:pos0'.
    Called inside the pileup loop — must be fast.
    """
    flags = 0
    in_dbsnp  = key in dbsnp
    in_sample = key in sample_snps
    in_het    = key in het_positions
    n = in_dbsnp + in_sample + in_het
    if n == 0:          return PASS
    if n > 1:           return SNP_MULTI
    if in_dbsnp:        return SNP_KNOWN
    if in_sample:       return SNP_SAMPLE
    return SNP_HET


# ════════════════════════════════════════════════════════════════════════════
# Context helpers (unchanged from original taps_calling_fast.py)
# ════════════════════════════════════════════════════════════════════════════

def get_context(seq, pos):
    l = seq[pos-1] if pos > 0          else "N"
    r = seq[pos+1] if pos+1 < len(seq) else "N"
    return l + seq[pos] + r


def context_match(seq, pos, is_reverse, context_filter):
    if context_filter == "ALL":
        return True
    if not is_reverse:
        downstream = seq[pos+1] if pos+1 < len(seq) else "N"
        if context_filter == "CpG":
            return downstream == "G"
        if context_filter == "CHG":
            return downstream != "G" and (pos+2 < len(seq) and seq[pos+2] == "G")
        if context_filter == "CHH":
            h1 = seq[pos+1] if pos+1 < len(seq) else "N"
            h2 = seq[pos+2] if pos+2 < len(seq) else "N"
            return h1 != "G" and h2 != "G"
    else:
        if context_filter == "CpG":
            return pos > 0 and seq[pos-1] == "C"
        if context_filter == "CHG":
            return pos >= 2 and seq[pos-1] != "C" and seq[pos-2] == "C"
        if context_filter == "CHH":
            h1 = seq[pos-1] if pos >= 1 else "N"
            h2 = seq[pos-2] if pos >= 2 else "N"
            return h1 != "C" and h2 != "C"
    return True


# ════════════════════════════════════════════════════════════════════════════
# Chromosome-level pileup worker
# SNP check happens HERE — before counts are accumulated
# ════════════════════════════════════════════════════════════════════════════

def process_chromosome(job):
    """
    Per-chromosome pileup worker.

    SNP filtering is applied at the position level BEFORE accumulating
    mod/unmod counts. Flagged positions are recorded but not counted —
    they are returned as zero-count entries with their flag so the
    caller can log how many were filtered per chromosome if needed.

    Returns list of dicts: one per PASS position meeting min_cov.
    """
    (bam_path, fasta_path, chrom,
     min_qual, min_mapq, min_cov, context_filter,
     dbsnp, sample_snps, het_positions) = job

    bam   = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    try:
        ref_seq = fasta.fetch(chrom).upper()
    except (KeyError, ValueError):
        bam.close(); fasta.close()
        return []

    counts = defaultdict(lambda: {"mod": 0.0, "unmod": 0.0, "ctx": "", "snp_flag": PASS})
    snp_skipped = 0

    try:
        reads = bam.fetch(chrom)
    except (ValueError, KeyError):
        bam.close(); fasta.close()
        return []

    for read in reads:
        if read.is_unmapped or read.query_sequence is None:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        # Bowtie1 uses XA:i:N; standard SAM uses NH:i:N
        try:
            nh = read.get_tag("NH")
        except KeyError:
            try:
                nh = read.get_tag("XA")
            except KeyError:
                nh = 1
        weight = 1.0 / max(nh, 1)

        quals   = read.query_qualities
        reverse = read.is_reverse

        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos >= len(ref_seq):
                continue
            if quals is not None and quals[qpos] < min_qual:
                continue

            fb = ref_seq[rpos].upper()

            # ── SNP CHECK: before any counting ───────────────────────────
            # Check at the position level. This is the correct place:
            # we reject the entire genomic position, not individual reads.
            # A SNP affects all reads at that locus — filtering per-read
            # would still accumulate counts and produce a biased rate.
            if fb == "C" or fb == "G":
                key      = f"{chrom}:{rpos}"
                snp_flag = assign_snp_flag(key, dbsnp, sample_snps, het_positions)
                if snp_flag != PASS:
                    # Record flag for this position but skip counting
                    if rpos not in counts:
                        counts[rpos]["snp_flag"] = snp_flag
                    continue   # ← skip all read-level accumulation
            # ── End SNP check ─────────────────────────────────────────────

            rb = read.query_sequence[qpos].upper()

            if not reverse:
                if fb != "C":
                    continue
                if not context_match(ref_seq, rpos, False, context_filter):
                    continue
                counts[rpos]["ctx"] = get_context(ref_seq, rpos)
                if rb == "C":   counts[rpos]["unmod"] += weight
                elif rb == "T": counts[rpos]["mod"]   += weight
            else:
                if fb != "G":
                    continue
                if not context_match(ref_seq, rpos, True, context_filter):
                    continue
                counts[rpos]["ctx"] = get_context(ref_seq, rpos)
                if rb == "G":   counts[rpos]["unmod"] += weight
                elif rb == "A": counts[rpos]["mod"]   += weight

    bam.close()
    fasta.close()

    rows = []
    for pos, v in counts.items():
        # Skip positions that were SNP-flagged (no counts accumulated)
        if v["snp_flag"] != PASS:
            continue
        total = v["mod"] + v["unmod"]
        if total < min_cov:
            continue
        rows.append({
            "chrom":       chrom,
            "start":       pos,
            "end":         pos + 1,
            "context":     v["ctx"],
            "mod_count":   round(v["mod"],   2),
            "unmod_count": round(v["unmod"], 2),
            "coverage":    round(total,      2),
            "mod_rate":    round(v["mod"] / total, 4),
            "snp_flag":    PASS,
        })
    return rows


# ════════════════════════════════════════════════════════════════════════════
# Statistical testing — runs ONLY on PASS sites
# ════════════════════════════════════════════════════════════════════════════

def binomial_test_bh(rows, background_rate):
    """
    One-sided binomial test per site (H1: mod_rate > background_rate).
    BH FDR correction applied across ALL tested sites simultaneously.

    Because SNP sites were removed before this step, the correction pool
    contains only genuine candidate modification sites. This is important:
    highly significant SNP positions (mod_rate ~0.5, very high coverage)
    would otherwise dominate the FDR correction and distort padj for
    real low-coverage modification sites.
    """
    if not rows:
        return rows

    pvalues = [
        stats.binomtest(
            int(round(r["mod_count"])),
            int(round(r["coverage"])),
            background_rate,
            alternative="greater",
        ).pvalue
        for r in rows
    ]

    _, padj, _, _ = multipletests(pvalues, method="fdr_bh")

    for r, pval, adj in zip(rows, pvalues, padj):
        r["pvalue"] = round(float(pval), 6)
        r["padj"]   = round(float(adj),  6)

    return rows


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()

    # Load SNP resources once — passed to every worker via job tuple
    dbsnp                    = load_dbsnp(args.dbsnp_vcf)
    sample_snps, het_positions = load_sample_snp_bed(args.sample_snp_bed, args.het_threshold)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        chroms = [sq["SN"] for sq in bam.header.to_dict()["SQ"]]

    print(f"[taps] BAM        : {args.bam}",        flush=True)
    print(f"[taps] Chroms     : {len(chroms)}",     flush=True)
    print(f"[taps] Threads    : {args.threads}",    flush=True)
    print(f"[taps] Context    : {args.context}",    flush=True)
    print(f"[taps] Min cov    : {args.min_cov}",    flush=True)
    print(f"[taps] Min mapq   : {args.min_mapq}",   flush=True)
    print(f"[taps] BG rate    : {args.background_rate}", flush=True)
    print(f"[taps] SNP filter : BEFORE counting (correct order)", flush=True)

    # Build jobs — SNP sets included so each worker can filter independently
    jobs = [
        (args.bam, args.fasta, chrom,
         args.min_qual, args.min_mapq, args.min_cov, args.context,
         dbsnp, sample_snps, het_positions)
        for chrom in chroms
    ]

    all_rows = []
    with Pool(processes=args.threads) as pool:
        for i, result in enumerate(pool.imap_unordered(process_chromosome, jobs), 1):
            all_rows.extend(result)
            if i % 10 == 0 or i == len(chroms):
                print(f"[taps] {i}/{len(chroms)} chroms  ({len(all_rows)} PASS sites so far)",
                      flush=True)

    print(f"[taps] PASS sites after SNP filter : {len(all_rows):,}", flush=True)

    # Binomial test + BH correction — clean pool, no SNPs
    all_rows = binomial_test_bh(all_rows, args.background_rate)

    cols = ["chrom", "start", "end", "context",
            "mod_count", "unmod_count", "coverage", "mod_rate",
            "pvalue", "padj", "snp_flag"]

    df = pd.DataFrame(all_rows, columns=cols)
    df.sort_values(["chrom", "start"], inplace=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[taps] Written : {len(df):,} sites → {args.out}", flush=True)


if __name__ == "__main__":
    main()
