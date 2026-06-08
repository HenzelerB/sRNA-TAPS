# -*- coding: utf-8 -*-
"""
05_annotate_biotype.py
Annotate BAM reads by RNA biotype using Ensembl GTF.
Splits BAM into per-biotype BAMs for separate TAPS m5C calling.

Priority (highest first): miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other

Improvements over original:
  - Interval-index lookup (binary search) instead of linear scan
    → ~10x faster on chromosomes with many gene features
  - Handles 'chr' prefix mismatch between BAM and GTF
  - Writes biotype summary TSV with reads + percent per biotype
"""

import pysam
import argparse
import bisect
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bam",     required=True,  help="Sorted indexed BAM")
    p.add_argument("--gtf",     required=True,  help="Ensembl GTF")
    p.add_argument("--out_dir", required=True,  help="Output directory")
    p.add_argument("--sample",  required=True,  help="Sample name")
    return p.parse_args()


BIOTYPE_PRIORITY = {
    "miRNA" : 1, "tRNA"  : 2, "piRNA" : 3,
    "snoRNA": 4, "snRNA" : 5, "rRNA"  : 6,
    "lncRNA": 7, "other" : 99,
}

BIOTYPE_MAP = {
    "miRNA": "miRNA", "pre_miRNA": "miRNA",
    "tRNA": "tRNA", "Mt_tRNA": "tRNA",
    "piRNA": "piRNA",
    "snoRNA": "snoRNA", "scaRNA": "snoRNA",
    "snRNA": "snRNA",
    "rRNA": "rRNA", "Mt_rRNA": "rRNA", "rRNA_pseudogene": "rRNA",
    "lncRNA": "lncRNA", "lincRNA": "lncRNA",
    "vault_RNA": "other", "Y_RNA": "other", "misc_RNA": "other",
}

BIOTYPE_DIRS = ["miRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "rRNA", "lncRNA", "other"]


def parse_gtf(gtf_path):
    """
    Parse GTF gene-level features into a dict:
        chrom -> sorted list of (start, end, biotype)
    Also builds a start-position index per chrom for binary search.
    Returns: intervals dict, starts dict (both keyed by chrom).
    """
    print(f"Parsing GTF: {gtf_path}", flush=True)
    raw = defaultdict(list)
    loaded = 0

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue

            chrom  = cols[0]
            start  = int(cols[3]) - 1   # 0-based
            end    = int(cols[4])
            attrs  = cols[8]

            biotype = "other"
            for tok in attrs.split(";"):
                tok = tok.strip()
                if tok.startswith("gene_biotype"):
                    biotype = tok.split('"')[1]
                    break

            mapped = BIOTYPE_MAP.get(biotype, "other")
            raw[chrom].append((start, end, mapped))
            loaded += 1

    # Sort and build start-index per chrom
    intervals = {}
    starts    = {}
    for chrom, feats in raw.items():
        feats.sort(key=lambda x: x[0])
        intervals[chrom] = feats
        starts[chrom]    = [f[0] for f in feats]

    print(f"  Loaded {loaded:,} gene features across {len(intervals)} chroms", flush=True)
    return intervals, starts


def find_biotype(chrom, pos, intervals, starts, chrom_aliases):
    """
    Binary-search for overlapping features at pos.
    Tries chrom aliases to handle 'chr1' vs '1' mismatches.
    Returns biotype string.
    """
    for c in chrom_aliases.get(chrom, [chrom]):
        if c not in intervals:
            continue
        feat_list  = intervals[c]
        start_list = starts[c]

        # Find rightmost feature that starts at or before pos
        idx = bisect.bisect_right(start_list, pos) - 1
        if idx < 0:
            continue

        best_priority = 999
        best_biotype  = "other"

        # Check backwards from idx while features could overlap pos
        i = idx
        while i >= 0 and start_list[i] >= pos - 200_000:
            feat_start, feat_end, biotype = feat_list[i]
            if feat_start <= pos <= feat_end:
                priority = BIOTYPE_PRIORITY.get(biotype, 99)
                if priority < best_priority:
                    best_priority = priority
                    best_biotype  = biotype
            i -= 1

        return best_biotype

    return "other"


def build_chrom_aliases(bam_chroms, gtf_chroms):
    """
    Map BAM chromosome names to GTF names.
    Handles 'chr1' ↔ '1' and 'chrM' ↔ 'MT' mismatches.
    """
    gtf_set = set(gtf_chroms)
    aliases = {}
    for c in bam_chroms:
        candidates = [c]
        if c.startswith("chr"):
            stripped = c[3:]
            candidates.append(stripped)
            if stripped == "M":
                candidates.append("MT")
        else:
            candidates.append("chr" + c)
        aliases[c] = [x for x in candidates if x in gtf_set]
    return aliases


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    intervals, starts = parse_gtf(args.gtf)

    bam_in = pysam.AlignmentFile(args.bam, "rb")
    header = bam_in.header

    bam_chroms    = [sq["SN"] for sq in header.to_dict()["SQ"]]
    chrom_aliases = build_chrom_aliases(bam_chroms, list(intervals.keys()))

    # Open per-biotype output BAMs
    out_bams = {}
    for bt in BIOTYPE_DIRS:
        bt_dir = out_dir / bt
        bt_dir.mkdir(parents=True, exist_ok=True)
        out_bams[bt] = pysam.AlignmentFile(
            str(bt_dir / f"{args.sample}_{bt}.bam"), "wb", header=header
        )

    counts = defaultdict(int)
    total  = 0

    print(f"Annotating: {args.bam}", flush=True)
    for read in bam_in.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        total += 1
        chrom   = read.reference_name
        pos     = read.reference_start
        biotype = find_biotype(chrom, pos, intervals, starts, chrom_aliases)
        counts[biotype] += 1
        out_bams[biotype].write(read)
        if total % 1_000_000 == 0:
            print(f"  {total:,} reads processed...", flush=True)

    bam_in.close()
    for bt, bam_out in out_bams.items():
        bam_out.close()

    # Sort, index, clean up
    print("Sorting and indexing biotype BAMs...", flush=True)
    for bt in BIOTYPE_DIRS:
        unsorted = str(out_dir / bt / f"{args.sample}_{bt}.bam")
        sorted_  = str(out_dir / bt / f"{args.sample}_{bt}.sorted.bam")
        if counts.get(bt, 0) > 0:
            pysam.sort("-@", "4", "-o", sorted_, unsorted)
            pysam.index(sorted_)
        else:
            # Create empty sorted BAM so Snakemake finds expected output
            bam_in = pysam.AlignmentFile(args.bam, "rb")
            empty  = pysam.AlignmentFile(sorted_, "wb", header=bam_in.header)
            empty.close()
            bam_in.close()
            pysam.index(sorted_)
        Path(unsorted).unlink(missing_ok=True)

    # Write summary
    summary_path = out_dir / f"{args.sample}_biotype_summary.txt"
    print(f"\n{'='*55}\nBiotype composition: {args.sample}\n{'='*55}")
    print(f"  {'Biotype':<12} {'Reads':>12} {'Percent':>9}")
    print(f"  {'-'*35}")
    with open(summary_path, "w") as sf:
        sf.write("biotype\treads\tpercent\n")
        for bt in sorted(counts, key=lambda x: counts[x], reverse=True):
            pct = counts[bt] / total * 100 if total > 0 else 0
            print(f"  {bt:<12} {counts[bt]:>12,} {pct:>8.1f}%")
            sf.write(f"{bt}\t{counts[bt]}\t{pct:.2f}\n")
    print(f"  {'TOTAL':<12} {total:>12,}   100.0%")
    print(f"\n  Summary: {summary_path}", flush=True)


if __name__ == "__main__":
    main()
