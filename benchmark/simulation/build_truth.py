#!/usr/bin/env python3
"""
Stage 2b -- build the ground-truth m5C site table for the sRNA-TAPS test-set
rebuild (Fix A).

For each vetted locus:
  1. fetch its plus-strand genomic sequence from the indexed FASTA
  2. reconstruct the *transcript* sequence (revcomp for minus-strand loci, so it
     matches what the simulator must emit under bowtie `--norc`)
  3. choose a deterministic, gene_id-seeded subset of its C positions to plant
     as m5C ground-truth sites
  4. map each planted C back to a strand-aware GENOMIC coordinate
  5. assign it a treat-condition planted rate (default 0.05-0.95)

Emits truth.tsv with columns:
    chrom  genomic_pos  strand  biotype  gene_id  planted_rate

Coordinate contract (matches align step `-v 2 --norc -k 10 --best --strata`):
  * vetted TSV start1/end1 are 1-based, inclusive.
  * '+' locus: transcript C at offset i  -> genomic_pos = start1 + i, ref base C.
  * '-' locus: transcript C at offset i  -> genomic_pos = end1 - i,  ref base G
    (a C on the minus strand is a G in the plus-strand reference).

This script is READ-ONLY against the genome and the probe lists; it only writes
truth.tsv (and its parent dir). It does not touch sRNA_TAPS_test/.
"""
import argparse
import csv
import hashlib
import os
import random
import sys

import pysam

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def stable_rng(gene_id: str, global_seed: int) -> random.Random:
    """Deterministic per-locus RNG so reruns are byte-identical."""
    h = hashlib.sha256(f"{global_seed}:{gene_id}".encode()).hexdigest()
    return random.Random(int(h[:16], 16))


def load_loci(path: str):
    """Yield dicts for each locus. Tolerates an optional header line and either
    tab- or whitespace-delimited rows. Expected fields, in order:
        chrom  start1  end1  strand  gene_id  [n_best_hits]
    """
    with open(path) as fh:
        first = fh.readline()
        fh.seek(0)
        # header if the 2nd column isn't an integer
        has_header = True
        parts = first.rstrip("\n").split("\t") if "\t" in first else first.split()
        if len(parts) >= 3:
            try:
                int(parts[1]); int(parts[2])
                has_header = False
            except ValueError:
                has_header = True
        reader = (l.rstrip("\n") for l in fh)
        if has_header:
            next(reader, None)
        for ln, line in enumerate(reader, start=2 if has_header else 1):
            if not line.strip():
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 5:
                raise ValueError(f"{path}:{ln}: expected >=5 cols, got {len(f)}: {line!r}")
            yield {
                "chrom": f[0],
                "start1": int(f[1]),
                "end1": int(f[2]),
                "strand": f[3],
                "gene_id": f[4],
                "n_best_hits": int(f[5]) if len(f) > 5 and f[5].lstrip("-").isdigit() else None,
            }


def plant_sites_for_locus(fa, locus, biotype, frac, min_sites, max_sites,
                          treat_min, treat_max, global_seed):
    chrom = locus["chrom"]; start1 = locus["start1"]; end1 = locus["end1"]
    strand = locus["strand"]; gene_id = locus["gene_id"]

    if end1 < start1:
        raise ValueError(f"{gene_id}: end1 < start1 ({start1}..{end1})")

    # pysam fetch is 0-based, half-open. start1-1 .. end1 gives the inclusive span.
    plus = fa.fetch(chrom, start1 - 1, end1).upper()
    L = end1 - start1 + 1
    if len(plus) != L:
        raise ValueError(
            f"{gene_id}: fetched {len(plus)} bp, expected {L} "
            f"({chrom}:{start1}-{end1}) -- coord/index mismatch")

    transcript = plus if strand == "+" else revcomp(plus)

    # C positions on the transcript (5'->3' as the read is emitted)
    c_offsets = [i for i, b in enumerate(transcript) if b == "C"]
    if not c_offsets:
        return [], 0  # no C -> nothing to plant

    rng = stable_rng(gene_id, global_seed)
    k = round(frac * len(c_offsets))
    k = max(min_sites, min(k, max_sites, len(c_offsets)))
    chosen = sorted(rng.sample(c_offsets, k))

    rows = []
    for i in chosen:
        if strand == "+":
            gpos = start1 + i
            ref = plus[i]
            assert ref == "C", f"{gene_id}: ref base at +{gpos} is {ref}, not C"
        else:
            gpos = end1 - i
            ref = plus[L - 1 - i]
            assert ref == "G", f"{gene_id}: ref base at -{gpos} is {ref}, not G"
        rate = round(rng.uniform(treat_min, treat_max), 4)
        rows.append((chrom, gpos, strand, biotype, gene_id, rate))
    return rows, len(c_offsets)


def main():
    ap = argparse.ArgumentParser(description="Stage 2b: build truth.tsv for sRNA-TAPS test rebuild")
    ap.add_argument("--probe-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "probes"),
                    help="dir holding probe_<biotype>_vetted.tsv")
    ap.add_argument("--fasta", required=True,
                    help="indexed genome FASTA (.fai must exist alongside)")
    ap.add_argument("--biotypes", nargs="+",
                    default=["miRNA", "snoRNA", "snRNA", "rRNA", "tRNA"])
    ap.add_argument("--out", default="truth_shared.tsv")
    ap.add_argument("--frac", type=float, default=0.25,
                    help="target fraction of a locus's Cs to plant as m5C")
    ap.add_argument("--min-sites", type=int, default=1)
    ap.add_argument("--max-sites", type=int, default=3)
    ap.add_argument("--treat-min", type=float, default=0.05)
    ap.add_argument("--treat-max", type=float, default=0.95)
    ap.add_argument("--seed", type=int, default=20260625)
    args = ap.parse_args()

    if not os.path.exists(args.fasta + ".fai"):
        sys.exit(f"ERROR: FASTA index not found: {args.fasta}.fai  (run: samtools faidx {args.fasta})")

    fa = pysam.FastaFile(args.fasta)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    all_rows = []
    print(f"{'biotype':<8} {'loci':>5} {'with_site':>9} {'no_C':>5} {'sites':>6}", file=sys.stderr)
    print("-" * 40, file=sys.stderr)
    for biotype in args.biotypes:
        probe = os.path.join(args.probe_dir, f"probe_{biotype}_vetted.tsv")
        if not os.path.exists(probe):
            sys.exit(f"ERROR: missing probe list: {probe}")
        n_loci = n_with = n_noC = n_sites = 0
        for locus in load_loci(probe):
            n_loci += 1
            rows, n_c = plant_sites_for_locus(
                fa, locus, biotype, args.frac, args.min_sites, args.max_sites,
                args.treat_min, args.treat_max, args.seed)
            if n_c == 0:
                n_noC += 1
            if rows:
                n_with += 1
                n_sites += len(rows)
                all_rows.extend(rows)
        print(f"{biotype:<8} {n_loci:>5} {n_with:>9} {n_noC:>5} {n_sites:>6}", file=sys.stderr)

    all_rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["chrom", "genomic_pos", "strand", "biotype", "gene_id", "planted_rate"])
        w.writerows(all_rows)

    print("-" * 40, file=sys.stderr)
    print(f"wrote {len(all_rows)} planted m5C sites -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
