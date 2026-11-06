#!/usr/bin/env python3
"""
simulate_taps_srna_v2.py
========================
Stage 2c -- generate synthetic TAPS small-RNA FASTQs for the sRNA-TAPS test-set
rebuild (Fix A).

Unlike v1 (which emitted hand-typed RNA templates that landed on the wrong
genomic loci and were mis-assigned by the biotype splitter), this version draws
reads from the ACTUAL genomic transcript of each vetted locus and plants m5C at
exactly the positions recorded in truth.tsv (Stage 2b). Because every read is a
substring of a correctly-annotated locus, the biotype splitter assigns it to the
intended biotype.

Produces 100 samples (HEK only), identical naming / samples.tsv schema to v1:
    no-treat_Ctrl_HEK_R1..R33   background C->T only
    pb_Ctrl_HEK_R1..R33         background + slight PB signal at m5C sites
    treat_HEK_R1..R34           genuine TAPS signal at m5C sites

Key contracts:
  * truth.tsv is the SINGLE source of truth: which C positions are m5C and, for
    the treat condition, the planted rate. The simulator consumes it; it does
    not re-derive truth.
  * Strand / coordinate math is IMPORTED from build_truth.py (revcomp, load_loci)
    so the two stages cannot drift. Minus-strand loci are emitted as the revcomp
    transcript, honoring bowtie `--norc` (a read that is a substring of the
    transcript is what maps back to the vetted locus).
  * Chemistry RANGES (updated for the full-spectrum truth set, planted_rate
    now spanning U(0.05, 0.95) rather than the original U(0.60, 0.90) band):
        no_treat       : every C  ~ U(0.01, 0.03)
        pb_ctrl        : "shared" m5C ~ planted_rate x U(0.05, 0.15) --
                         scales with the site's own signal instead of the old
                         fixed U(0.60, 0.70), which would otherwise sit ABOVE
                         a widened low treat rate and invert the intended
                         treat > pb_ctrl > no_treat ordering.
                         "discovery" m5C and all non-m5C C's ~ U(0.01, 0.03)
        treat          : m5C       = planted_rate (truth.tsv); other C ~ U(0.01,0.03)
    With --redraw-treat-rate, treat m5C reverts to per-sample U(0.60,0.90) (v1).
  * Per-replicate biological variance: each of the 33/34 replicates of a
    condition now draws its OWN rate per site from Beta(mu*kappa,
    (1-mu)*kappa) around the base rate above (mu), instead of all replicates
    sharing one fixed rate with only binomial (per-read) sampling noise on
    top -- the old behaviour understated real replicate-to-replicate
    dispersion (verified: cross-replicate SD matched pure sqrt(p(1-p)/n)
    almost exactly, i.e. zero biological variance). Default --replicate-kappa
    50 targets ~12% CV, in line with published bisulfite/TAPS biological
    replicate variability; pass a negative value to disable and restore the
    old single-shared-rate behaviour.
  * Read allocation is proportional to locus length -> uniform per-position depth,
    so every planted site is equally testable across biotypes (fixes v1 miRNA
    dominance + small-biotype starvation).

Usage:
    python3 simulate_taps_srna_v2.py \
        --probe-dir sRNA_TAPS_test \
        --fasta     sRNA_TAPS_bio/04a.genome/Homo_sapiens.GRCh38.dna.toplevel.fa \
        --truth     sRNA_TAPS_test_v2/truth.tsv \
        --outdir    sRNA_TAPS_test_v2/fastq \
        --reads     1000000
"""

import argparse
import gzip
import os
import random
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# Import the strand/coordinate contract from Stage 2b so the two cannot diverge.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_truth import revcomp, load_loci  # noqa: E402

ADAPTER_SEQ = "TGGAATTCTCGGGTGCCAAGG"   # TruSeq small RNA 3' adapter (unchanged)
BASE_ERROR_RATE = 0.005
BASES = "ACGT"

DEFAULT_BIOTYPES = ["miRNA", "snoRNA", "snRNA", "rRNA", "tRNA"]
# per-biotype read-length windows (min, max); capped at locus length per read
READ_LEN_RANGES = {
    "miRNA":  (18, 24),
    "snoRNA": (18, 28),
    "snRNA":  (18, 28),
    "rRNA":   (18, 26),
    "tRNA":   (18, 28),
}


# ── low-level read mechanics (unchanged from v1) ──────────────────────────────
def mutate_base(base, rng, error_rate):
    if rng.random() < error_rate:
        return rng.choice([b for b in BASES if b != base])
    return base


def add_sequencing_errors(seq, rng, error_rate=BASE_ERROR_RATE):
    return "".join(mutate_base(b, rng, error_rate) for b in seq)


def generate_quality_string(length, rng, min_q=30, max_q=40):
    quals = [rng.randint(min_q, max_q) for _ in range(length)]
    for i in range(min(3, length)):
        quals[i] = max(20, quals[i] - rng.randint(5, 10))
    for i in range(max(0, length - 3), length):
        quals[i] = max(20, quals[i] - rng.randint(5, 10))
    return "".join(chr(q + 33) for q in quals)


# ── reference construction from vetted loci + truth.tsv ───────────────────────
def load_truth(path):
    """gene_id -> list of (genomic_pos:int, strand, planted_rate:float)."""
    by_gene = defaultdict(list)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for need in ("genomic_pos", "strand", "gene_id", "planted_rate"):
            if need not in idx:
                sys.exit(f"ERROR: truth.tsv missing column '{need}' (header: {header})")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            by_gene[f[idx["gene_id"]]].append((
                int(f[idx["genomic_pos"]]),
                f[idx["strand"]],
                float(f[idx["planted_rate"]]),
                f[idx["truth_class"]] if "truth_class" in idx else "shared",
            ))
    return by_gene


def build_reference(probe_dir, fasta, biotypes, truth_path):
    """Return dict biotype -> list of locus dicts:
        {gene_id, chrom, start1, end1, strand, transcript(DNA str),
         length, m5c: {offset: planted_rate}}
    """
    if not os.path.exists(fasta + ".fai"):
        sys.exit(f"ERROR: FASTA index not found: {fasta}.fai")
    fa = pysam.FastaFile(fasta)
    truth = load_truth(truth_path)

    reference = {}
    seen_genes = {}
    for biotype in biotypes:
        probe = os.path.join(probe_dir, f"probe_{biotype}_vetted.tsv")
        if not os.path.exists(probe):
            sys.exit(f"ERROR: missing probe list: {probe}")
        loci = []
        for locus in load_loci(probe):
            gid = locus["gene_id"]
            if gid in seen_genes:
                sys.exit(f"ERROR: gene_id {gid} appears in both "
                         f"{seen_genes[gid]} and {biotype} probe lists")
            seen_genes[gid] = biotype

            chrom, s1, e1, strand = locus["chrom"], locus["start1"], locus["end1"], locus["strand"]
            plus = fa.fetch(chrom, s1 - 1, e1).upper()
            L = e1 - s1 + 1
            if len(plus) != L:
                sys.exit(f"ERROR: {gid}: fetched {len(plus)} bp != {L} ({chrom}:{s1}-{e1})")
            transcript = plus if strand == "+" else revcomp(plus)

            m5c = {}
            for gpos, tstrand, rate, truth_class in truth.get(gid, []):
                off = (gpos - s1) if strand == "+" else (e1 - gpos)
                if not (0 <= off < L):
                    sys.exit(f"ERROR: {gid}: truth pos {gpos} -> offset {off} out of range [0,{L})")
                if transcript[off] != "C":
                    sys.exit(f"ERROR: {gid}: truth pos {gpos} -> transcript[{off}]="
                             f"{transcript[off]}, expected C (strand mismatch with 2b)")
                m5c[off] = (rate, truth_class)

            loci.append({
                "gene_id": gid, "chrom": chrom, "start1": s1, "end1": e1,
                "strand": strand, "transcript": transcript, "length": L, "m5c": m5c,
            })
        reference[biotype] = loci
    return reference


def precompute_site_rates(reference, master_seed, redraw_treat):
    """Per (locus, condition), fix a BASE (mean) rate for every C offset.

    This used to be the rate every one of the 33/34 replicates simulated
    from directly, which meant the only variation between replicates was
    binomial sampling noise from finite read depth (empirically verified:
    observed cross-replicate SD matched sqrt(p(1-p)/n) almost exactly) --
    unlike real biological replicates, which disperse around the true rate
    (batch effects, extraction/conversion efficiency, day-to-day drift) on
    top of that sampling noise. This base rate is now only the MEAN each
    replicate's own rate is drawn around -- see draw_replicate_rates().

    Stored on each locus as locus['base_rates'][condition] = {offset: rate}.
    """
    for biotype, loci in reference.items():
        for locus in loci:
            c_offsets = [i for i, b in enumerate(locus["transcript"]) if b == "C"]
            m5c = locus["m5c"]
            rates = {}
            for cond in ("no_treat", "pb_ctrl", "treat"):
                rng = random.Random(f"{master_seed}:{locus['gene_id']}:{cond}")
                d = {}
                for o in c_offsets:
                    is_m5c = o in m5c
                    truth_class = m5c[o][1] if is_m5c else ""
                    if cond == "no_treat":
                        d[o] = rng.uniform(0.01, 0.03)
                    elif cond == "pb_ctrl":
                        # Scale with the site's own treat rate rather than a
                        # fixed band -- a fixed U(0.60,0.70) would sit ABOVE
                        # low planted rates once the truth set spans down to
                        # 0.05, inverting treat > pb_ctrl > no_treat.
                        d[o] = m5c[o][0] * rng.uniform(0.05, 0.15) if truth_class == "shared" else rng.uniform(0.01, 0.03)
                    else:  # treat
                        if is_m5c:
                            d[o] = rng.uniform(0.60, 0.90) if redraw_treat else m5c[o][0]
                        else:
                            d[o] = rng.uniform(0.01, 0.03)
                rates[cond] = d
            locus["base_rates"] = rates


def draw_replicate_rates(locus, condition, sample_name, master_seed, kappa):
    """This replicate's own per-offset rate, drawn from
    Beta(mu*kappa, (1-mu)*kappa) around the locus's base (mean) rate --
    standard beta-binomial overdispersion so replicates vary biologically,
    not just from per-read sampling noise. Deterministic per
    (locus, condition, sample) so reruns are reproducible.
    kappa=None disables this (returns the base rates unchanged).
    """
    base = locus["base_rates"][condition]
    if kappa is None:
        return base
    rng = random.Random(f"{master_seed}:{locus['gene_id']}:{condition}:{sample_name}:rep")
    eps = 1e-3
    out = {}
    for o, mu in base.items():
        mu_c = min(max(mu, eps), 1 - eps)
        out[o] = rng.betavariate(mu_c * kappa, (1 - mu_c) * kappa)
    return out


# ── read generation ───────────────────────────────────────────────────────────
def generate_read(locus, site_rates, rng, read_len_range):
    t = locus["transcript"]
    L = locus["length"]
    lo, hi = read_len_range
    read_len = rng.randint(lo, min(hi, L)) if L >= lo else L
    start = rng.randint(0, L - read_len)

    out = []
    for j in range(read_len):
        o = start + j
        b = t[o]
        if b == "C" and rng.random() < site_rates.get(o, 0.02):
            b = "T"            # TAPS: modified/converted C reads as T
        out.append(b)
    subseq = add_sequencing_errors("".join(out), rng)

    if rng.random() < 0.5:                       # adapter readthrough on ~half
        alen = rng.randint(4, len(ADAPTER_SEQ))
        subseq = subseq + ADAPTER_SEQ[:alen]
    subseq = subseq[:50]
    qual = generate_quality_string(len(subseq), rng)
    return subseq, qual


def allocate_reads(reference, total_reads, min_reads_per_locus):
    """Reads per locus proportional to locus length -> uniform per-position depth.
    Returns {gene_id: n_reads}. Guarantees >= min_reads_per_locus each.
    """
    lengths = {l["gene_id"]: l["length"]
               for loci in reference.values() for l in loci}
    total_len = sum(lengths.values())
    alloc = {}
    for gid, L in lengths.items():
        alloc[gid] = max(min_reads_per_locus, round(total_reads * L / total_len))
    return alloc


def generate_sample_fastq(sample_name, condition, alloc, reference, outdir, sample_rng,
                          master_seed, replicate_kappa):
    outfile = outdir / f"{sample_name}.fastq.gz"
    counter = 0
    with gzip.open(outfile, "wt") as fout:
        for biotype, loci in reference.items():
            rng_range = READ_LEN_RANGES.get(biotype, (18, 26))
            for locus in loci:
                site_rates = draw_replicate_rates(locus, condition, sample_name, master_seed, replicate_kappa)
                n = alloc[locus["gene_id"]]
                for _ in range(n):
                    seq, qual = generate_read(locus, site_rates, sample_rng, rng_range)
                    rid = f"@{sample_name}_{biotype}_{locus['gene_id']}_{counter}"
                    fout.write(f"{rid}\n{seq}\n+\n{qual}\n")
                    counter += 1
    return outfile, counter


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="Stage 2c: simulate TAPS sRNA reads from vetted loci")
    p.add_argument("--probe-dir", required=True, help="dir with probe_<biotype>_vetted.tsv")
    p.add_argument("--fasta", required=True, help="indexed genome FASTA (.fai alongside)")
    p.add_argument("--truth", required=True, help="truth.tsv from Stage 2b")
    p.add_argument("--outdir", required=True)
    p.add_argument("--reads", type=int, default=1_000_000, help="approx reads per sample")
    p.add_argument("--min-reads-per-locus", type=int, default=50)
    p.add_argument("--biotypes", nargs="+", default=DEFAULT_BIOTYPES)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--redraw-treat-rate", action="store_true",
                   help="v1 behaviour: redraw treat m5C rate U(0.60,0.90) instead of truth.tsv")
    p.add_argument("--replicate-kappa", type=float, default=50.0,
                   help="beta-binomial dispersion for per-replicate rate around each site's "
                        "base rate (~12%% CV at kappa=50); pass a negative value to disable "
                        "(replicates then share one fixed rate, only binomial sampling noise)")
    args = p.parse_args()
    replicate_kappa = None if args.replicate_kappa < 0 else args.replicate_kappa

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reference = build_reference(args.probe_dir, args.fasta, args.biotypes, args.truth)
    precompute_site_rates(reference, args.seed, args.redraw_treat_rate)
    alloc = allocate_reads(reference, args.reads, args.min_reads_per_locus)

    # per-biotype depth report (verify BEFORE the pipeline run)
    print(f"{'biotype':<8} {'loci':>5} {'bp':>7} {'reads/sample':>13} {'~depth/pos':>11}", file=sys.stderr)
    print("-" * 50, file=sys.stderr)
    grand_reads = 0
    for biotype, loci in reference.items():
        bp = sum(l["length"] for l in loci)
        reads = sum(alloc[l["gene_id"]] for l in loci)
        grand_reads += reads
        avg_rl = sum(READ_LEN_RANGES.get(biotype, (18, 26))) / 2
        depth = reads * avg_rl / bp if bp else 0
        print(f"{biotype:<8} {len(loci):>5} {bp:>7} {reads:>13,} {depth:>11.0f}", file=sys.stderr)
    print("-" * 50, file=sys.stderr)
    print(f"~{grand_reads:,} reads/sample x 100 samples", file=sys.stderr)

    samples = []
    for condition, n_reps in [("no_treat", 33), ("pb_ctrl", 33), ("treat", 34)]:
        for r in range(1, n_reps + 1):
            if condition == "no_treat":
                name = f"no-treat_Ctrl_HEK_R{r}"
            elif condition == "pb_ctrl":
                name = f"pb_Ctrl_HEK_R{r}"
            else:
                name = f"treat_HEK_R{r}"
            samples.append((name, condition, r))

    master = random.Random(args.seed)
    print(f"\nWriting {len(samples)} samples -> {outdir}/", file=sys.stderr)
    for name, condition, rep in samples:
        sample_rng = random.Random(master.random())   # reproducible per-sample stream
        outfile, n = generate_sample_fastq(name, condition, alloc, reference, outdir, sample_rng,
                                           args.seed, replicate_kappa)
        size_mb = outfile.stat().st_size / 1e6
        print(f"  {name:<22} {n:>9,} reads  ({size_mb:5.1f} MB)", file=sys.stderr)

    tsv = outdir / "samples.tsv"
    with open(tsv, "w") as f:
        f.write("sample\tcondition\tcell_line\treplicate\tfastq\n")
        for name, condition, rep in samples:
            f.write(f"{name}\t{condition}\tHEK\tR{rep}\t{outdir / (name + '.fastq.gz')}\n")
    print(f"samples.tsv -> {tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
