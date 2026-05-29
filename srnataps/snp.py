# -*- coding: utf-8 -*-
"""
06_build_snp_blacklist.py
Build a per-cell-line SNP blacklist from no-treat BAMs.

The no-treat condition has no TAPS chemistry. Any C→T (plus strand) or
G→A (minus strand) at allele frequency >= min_af in these BAMs is a
germline variant in this cell line — not a modification signal.

HEK293 and Caco2 are immortalised cell lines with their own distinct
mutation profiles, many of which are not captured in population-derived
dbSNP. This script generates a cell-line-specific blacklist that is used
alongside dbSNP in 07_taps_calling.py (Layer 2 filtering).

Output: BED file with columns:
    chrom, start, end, ref, alt, allele_freq, coverage, cell_line

Usage:
    python 06_build_snp_blacklist.py \
        --bam no_treat_merged.bam \
        --fasta hg38.fa \
        --out sample_snps_HEK.bed \
        --min-af 0.20 \
        --min-cov 10 \
        --cell-line HEK
"""

import argparse
import sys
from pathlib import Path
import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam",        required=True)
    p.add_argument("--fasta",      required=True)
    p.add_argument("--out",        required=True)
    p.add_argument("--min-af",     type=float, default=0.20,
                   help="Minimum allele frequency to call a SNP (default: 0.20)")
    p.add_argument("--min-cov",    type=int,   default=10,
                   help="Minimum coverage at site (default: 10)")
    p.add_argument("--cell-line",  default="unknown")
    return p.parse_args()


def main():
    args = parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    fasta = pysam.FastaFile(args.fasta)
    bam   = pysam.AlignmentFile(args.bam, "rb")

    n_snps = 0
    print(f"[snp_blacklist] BAM        : {args.bam}", flush=True)
    print(f"[snp_blacklist] Cell line  : {args.cell_line}", flush=True)
    print(f"[snp_blacklist] Min AF     : {args.min_af}", flush=True)
    print(f"[snp_blacklist] Min cov    : {args.min_cov}", flush=True)

    with open(args.out, "w") as out:
        out.write("# Sample-specific C→T and G→A SNPs from no-treat BAM\n")
        out.write(f"# cell_line={args.cell_line}  min_af={args.min_af}  min_cov={args.min_cov}\n")
        out.write("chrom\tstart\tend\tref\talt\tallele_freq\tcoverage\tcell_line\n")

        for chrom in bam.references:
            try:
                seq = fasta.fetch(chrom).upper()
            except (KeyError, ValueError):
                continue

            for col in bam.pileup(
                chrom,
                min_mapping_quality=10,
                min_base_quality=20,
                ignore_overlaps=True,
                ignore_orphans=True,
            ):
                pos = col.reference_pos
                if pos >= len(seq):
                    continue
                ref = seq[pos]

                # Plus strand: C→T
                if ref == "C":
                    c = t = 0
                    for r in col.pileups:
                        if r.is_del or r.is_refskip:
                            continue
                        b = r.alignment.query_sequence[r.query_position].upper()
                        if b == "C": c += 1
                        elif b == "T": t += 1
                    cov = c + t
                    if cov >= args.min_cov and t / cov >= args.min_af:
                        af = t / cov
                        out.write(f"{chrom}\t{pos}\t{pos+1}\tC\tT\t{af:.4f}\t{cov}\t{args.cell_line}\n")
                        n_snps += 1

                # Minus strand: G→A (complement of C→T)
                elif ref == "G":
                    g = a = 0
                    for r in col.pileups:
                        if r.is_del or r.is_refskip:
                            continue
                        b = r.alignment.query_sequence[r.query_position].upper()
                        if b == "G": g += 1
                        elif b == "A": a += 1
                    cov = g + a
                    if cov >= args.min_cov and a / cov >= args.min_af:
                        af = a / cov
                        out.write(f"{chrom}\t{pos}\t{pos+1}\tG\tA\t{af:.4f}\t{cov}\t{args.cell_line}\n")
                        n_snps += 1

    bam.close()
    fasta.close()

    print(f"[snp_blacklist] SNPs found : {n_snps:,}", flush=True)
    print(f"[snp_blacklist] Output     : {args.out}", flush=True)


if __name__ == "__main__":
    main()
