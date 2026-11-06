#!/usr/bin/env python3
"""Add 100 discovery sites to the shared truth table to give the 738-site benchmark truth set.

Discovery sites are drawn from unplanted C positions of the same vetted loci (seed 20260827) with
planted rates from U(0.05, 0.95). Shared sites keep their rates from build_truth.py.
"""
import argparse
import csv
import os
import random

import pysam

from build_truth import load_loci, revcomp

BIOTYPES = ["miRNA", "snoRNA", "snRNA", "rRNA", "tRNA"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shared", default="truth_shared.tsv", help="output of build_truth.py")
    ap.add_argument("--probe-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "probes"))
    ap.add_argument("--fasta", required=True, help="indexed Ensembl GRCh38 dna.toplevel FASTA")
    ap.add_argument("--out", default="fig2c_truth_set.tsv")
    ap.add_argument("--n-discovery", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260827)
    args = ap.parse_args()

    with open(args.shared) as fh:
        shared = list(csv.DictReader(fh, delimiter="\t"))
    used = {(r["chrom"], r["genomic_pos"], r["strand"]) for r in shared}
    fa = pysam.FastaFile(args.fasta)
    candidates = []
    for biotype in BIOTYPES:
        for locus in load_loci(os.path.join(args.probe_dir, f"probe_{biotype}_vetted.tsv")):
            plus = fa.fetch(locus["chrom"], locus["start1"] - 1, locus["end1"]).upper()
            transcript = plus if locus["strand"] == "+" else revcomp(plus)
            for i, base in enumerate(transcript):
                if base != "C":
                    continue
                pos = locus["start1"] + i if locus["strand"] == "+" else locus["end1"] - i
                key = (locus["chrom"], str(pos), locus["strand"])
                if key not in used:
                    candidates.append((*key, biotype, locus["gene_id"]))

    rng = random.Random(args.seed)
    rng.shuffle(candidates)
    extra = candidates[: args.n_discovery]
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\r\n")
        writer.writerow(["chrom", "genomic_pos", "strand", "biotype", "gene_id", "planted_rate", "truth_class"])
        for r in shared:
            writer.writerow([r["chrom"], r["genomic_pos"], r["strand"], r["biotype"], r["gene_id"], r["planted_rate"], "shared"])
        for chrom, pos, strand, biotype, gene_id in extra:
            writer.writerow([chrom, pos, strand, biotype, gene_id, round(rng.uniform(0.05, 0.95), 4), "discovery"])
    print(f"wrote {len(shared)} shared + {len(extra)} discovery sites -> {args.out}")


if __name__ == "__main__":
    main()
