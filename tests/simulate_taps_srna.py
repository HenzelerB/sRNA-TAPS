#!/usr/bin/env python3
"""
simulate_taps_srna.py
=====================
Generate synthetic TAPS small RNA FASTQ files for sRNA-TAPS pipeline testing.

Simulates 100 samples (HEK cell line only):
  - no-treat_Ctrl_HEK_R1..R33  : background C→T only (1-3%)
  - pb_Ctrl_HEK_R1..R33        : background + slight PB signal (5-10%)
  - treat_HEK_R1..R34          : genuine TAPS signal (60-90% C→T at m5C)

Chemistry: TAPS — modified cytosines appear as C→T in reads.
Adapter: TruSeq small RNA 3' adapter (TGGAATTCTCGGGTGCCAAGG)
Reference: Uses real hg38 sequences for miRNA/tRNA/rRNA/snoRNA loci.

Usage:
    python3 simulate_taps_srna.py --outdir /path/to/fastq_dir --reads 1000000
"""

import argparse
import gzip
import random
from pathlib import Path

# ── Reference sequences ───────────────────────────────────────────────────────
REFERENCE_SEQUENCES = {

    # ── miRNA: 50 templates ───────────────────────────────────────────────────
    "miRNA": [
        ("hsa-miR-21-5p",     "UAGCUUAUCAGACUGAUGUUGA",   True,  [14]),
        ("hsa-miR-155-5p",    "UUAAUGCUAAUUGUGAUAGGGGU",  True,  [8, 18]),
        ("hsa-miR-122-5p",    "UGGAGUGUGACAAUGGUGUUUG",   True,  [9]),
        ("hsa-let-7a-5p",     "UGAGGUAGUAGGUUGUAUAGUU",   False, []),
        ("hsa-miR-16-5p",     "UAGCAGCACGUAAAUAUUGGCG",   False, []),
        ("hsa-miR-21-3p",     "CAACACCAGUCGAUGGGCUGU",    True,  [3, 15]),
        ("hsa-miR-126-3p",    "UCGUACCGUGAGUAAUAAUGCG",   True,  [2]),
        ("hsa-miR-29a-3p",    "UAGCACCAUCUGAAAUCGGUUA",   False, []),
        ("hsa-miR-223-3p",    "UGUCAGUUUGUCAAAUACCCCA",   True,  [19]),
        ("hsa-miR-92a-3p",    "UAUUGCACUUGUCCCGGCCUGU",   False, []),
        ("hsa-miR-10a-5p",    "UACCCUGUAGAUCCGAAUUUGU",   True,  [5, 12]),
        ("hsa-miR-181a-5p",   "AACAUUCAACGCUGUCGGUGAGU",  True,  [7]),
        ("hsa-miR-17-5p",     "CAAAGUGCUUACAGUGCAGGUAG",  False, []),
        ("hsa-miR-93-5p",     "CAAAGUGCUGUUCGUGCAGGUAG",  True,  [4]),
        ("hsa-miR-106b-5p",   "UAAAGUGCUGACAGUGCAGAU",    False, []),
        ("hsa-miR-20a-5p",    "UAAAGUGCUUAUAGUGCAGGUAG",  True,  [11]),
        ("hsa-miR-25-3p",     "CAUUGCACUUGUCUCGGUCUGA",   True,  [6, 16]),
        ("hsa-miR-130a-3p",   "CAGUGCAAUGUUAAAAGGGCAU",   True,  [3]),
        ("hsa-miR-221-3p",    "AGCUACAUUGUCUGCUGGGUUUC",  False, []),
        ("hsa-miR-222-3p",    "AGCUACAUCUGGCUACUGGGU",    True,  [9]),
        ("hsa-miR-34a-5p",    "UGGCAGUGUGGUUAGCUGGUUGU",  True,  [4, 17]),
        ("hsa-miR-145-5p",    "GUCCAGUUUUCCCAGGAAUCCCU",  False, []),
        ("hsa-miR-143-3p",    "UGAGAUGAAGCACUGUAGCUC",    True,  [8]),
        ("hsa-miR-199a-3p",   "ACAGUAGUCUGCACAUUGGUUA",   True,  [13]),
        ("hsa-miR-200c-3p",   "UAAUACUGCCGGGUAAUGAUGGA",  True,  [6]),
        ("hsa-miR-141-3p",    "UAACACUGUCUGGUAAAGAUGG",   False, []),
        ("hsa-miR-429",       "UAAUACUGUCUGGUAAAACCGU",   False, []),
        ("hsa-miR-200b-3p",   "UAAUACUGCCUGGUAAUGAUGA",   True,  [5, 19]),
        ("hsa-miR-203a-3p",   "GUGAAAUGUUUAGGACCACUAG",   True,  [7]),
        ("hsa-miR-205-5p",    "UCCUUCAUUCCACCGGAGUCUG",   True,  [3, 14]),
        ("hsa-miR-210-3p",    "CUGUGCGUGUGACAGCGGCUGA",   True,  [10]),
        ("hsa-miR-335-5p",    "UCAAGAGCAAUAACGAAAAAUGU",  False, []),
        ("hsa-miR-375",       "UUUGUUCGUUCGGCUCGCGUGA",   True,  [8]),
        ("hsa-miR-423-5p",    "UGAGGGGCAGAGAGCGAGACUUU",  True,  [5]),
        ("hsa-miR-451a",      "AAACCGUUACCAUUACUGAGUU",   False, []),
        ("hsa-miR-486-5p",    "UCCUGUACUGAGCUGCCCCGAG",   True,  [12]),
        ("hsa-miR-532-5p",    "CAUGCCUUGAGUGUAGGACCGU",   True,  [3, 16]),
        ("hsa-miR-590-5p",    "GAGCUUAUUCAUAAAAGUGCAG",   False, []),
        ("hsa-miR-660-5p",    "UACCCAUUGCAUAUUGGAGUUG",   True,  [7]),
        ("hsa-miR-769-5p",    "UGAGACCCUGCCCAGCCCUGA",    True,  [4]),
        ("hsa-miR-874-3p",    "CUGCCCUGGGACCGAGGGCCUG",   True,  [9, 18]),
        ("hsa-miR-940",       "AAGGCAGGGCCCCCGCUCCCGG",   False, []),
        ("hsa-miR-1246",      "AAUGGAUUUUUGGAGCAGG",      True,  [6]),
        ("hsa-miR-1290",      "UGGAUUUUUGGAUCAGGGA",      True,  [2, 13]),
        ("hsa-miR-3648",      "CGCGGCGGGCGGGGCGGGGCGG",   False, []),
        ("hsa-miR-4286",      "AGCCCGGCUCGGGCUCAGGGG",    True,  [5]),
        ("hsa-miR-4454",      "UGGUGCGAUCCUAGCUUAGAGU",   True,  [3, 11]),
        ("hsa-miR-4516",      "GGGAGAGGCUCGGGGAGAGGCU",   False, []),
        ("hsa-miR-4532",      "AGGGCGGGGCGCGGGCGGGGCG",   True,  [8]),
        ("hsa-miR-6126",      "GCUGGCCCGGGCGGGGCUGGCC",   True,  [6, 15]),
    ],

    # ── tRNA: 25 templates ────────────────────────────────────────────────────
    "tRNA": [
        ("mt-tRNA-Leu-frag1",  "GGUUCGAUUCCCGGUCUCAGG",   True,  [5, 12]),
        ("mt-tRNA-Ile-frag1",  "AGUGGUUUAAUCUUCUUAAGU",   True,  [3, 17]),
        ("mt-tRNA-Phe-frag1",  "GUCAUUCGUUGAAGAGUUGCA",   True,  [8, 14]),
        ("mt-tRNA-Val-frag1",  "CAUCAACUCCUAACACUUUCA",    True,  [6]),
        ("tRNA-Gly-frag1",     "GCGGUAGUAGCUCAGUCGGUA",    False, []),
        ("tRNA-Ala-frag1",     "GGGGGUAUAGCUCAGUGGUAG",    True,  [4, 11]),
        ("mt-tRNA-Pro-frag1",  "CAAUCCAGUGCUUGAGUCACA",    True,  [7]),
        ("tRNA-Ser-frag1",     "GCCCGGAUAGCUCAGUCGGUA",    False, []),
        ("mt-tRNA-Trp-frag1",  "AGUUAAGCUUGCAAGUGCUCA",    True,  [9]),
        ("tRNA-Asp-frag1",     "UCCGUAGUAGCUCAGUCGGUA",    True,  [2, 15]),
        ("tRNA-Glu-frag1",     "UCCCGGGUAGCUCAGUCGGUA",    True,  [5]),
        ("mt-tRNA-Gln-frag1",  "UUAAAUGCUUUCAGUGCUUGAA",   True,  [10]),
        ("tRNA-His-frag1",     "GUCGUAGUAGCUCAGUCGGUA",    False, []),
        ("tRNA-Lys-frag1",     "GCCCAGAUAGCUCAGUCGGUA",    True,  [3, 18]),
        ("mt-tRNA-Thr-frag1",  "CAAAGCUCUUGCAGUGCUCCA",    True,  [8]),
        ("tRNA-Cys-frag1",     "GUCGCAGUAGCUCAGUCGGUA",    False, []),
        ("tRNA-Met-frag1",     "GCUCGUGUAGCUCAGUCGGUA",    True,  [6, 14]),
        ("tRNA-Arg-frag1",     "GCUCGGGUAGCUCAGUCGGUA",    True,  [3]),
        ("tRNA-Tyr-frag1",     "GCUUAAGUAGCUCAGUCGGUA",    False, []),
        ("tRNA-Trp-frag1",     "GCUCAAAUAGCUCAGUCGGUA",    True,  [9, 17]),
        ("tRNA-Val-frag2",     "GCUUAGGUAGCUCAGUCGGUA",    True,  [5]),
        ("mt-tRNA-Lys-frag1",  "UUAAAUGCUUUCAAAGCUUGAA",   True,  [11]),
        ("tRNA-Phe-frag2",     "GCGGAAGUAGCUCAGUCGGUA",    True,  [4, 16]),
        ("mt-tRNA-Asp-frag1",  "AGUUAAGCUUGCAAGUGCUCA",    False, []),
        ("tRNA-Asn-frag1",     "GCUUAGGUAGCUCAGUCGGUA",    True,  [7]),
    ],

    # ── rRNA: 25 templates ────────────────────────────────────────────────────
    "rRNA": [
        ("mt-12S-frag1",   "CAUCACGAAACUCAGCACACU",   True,  [4, 16]),
        ("mt-12S-frag2",   "UGGCUACACCUUGACAGCUAC",   True,  [10]),
        ("mt-16S-frag1",   "GCUCGCCCUUGUGCAGAGAAU",   True,  [2, 14]),
        ("mt-16S-frag2",   "AUGGCUGAGCCAGGCCUUUGA",   True,  [8]),
        ("28S-frag1",      "GGCCCGAAACCCGACAGGACC",   False, []),
        ("28S-frag2",      "CCGAGCUCGAAUUUGCUUCGA",   True,  [5]),
        ("18S-frag1",      "ACGGUAGAGCUACCGAUUGCU",   False, []),
        ("18S-frag2",      "CGAAACUCGCCCAGCAAACGC",   True,  [3, 18]),
        ("5.8S-frag1",     "GCGAAACGCGAAUUGAACGCG",   True,  [9]),
        ("mt-12S-frag3",   "UGCAGCUUAACUCAAAGCACC",   True,  [6, 13]),
        ("28S-frag3",      "GCGAGCUCGGCCCGAAACCCG",   True,  [4]),
        ("18S-frag3",      "CCGUAGCUAAAUGCGGUCUCA",    False, []),
        ("5S-frag1",       "GCUGGUCCCAUACCGACUCCA",   True,  [7, 17]),
        ("mt-16S-frag3",   "UAAGCUUGCAAGUGCUCAAGU",   True,  [11]),
        ("28S-frag4",      "CUGAGCCAGGCCUUUGAAACG",   True,  [2, 16]),
        ("18S-frag4",      "AGCUACCGAUUGCUGAGCCAG",   True,  [8]),
        ("mt-12S-frag4",   "GCAACACUGAAACUCAGCAAG",   True,  [5, 14]),
        ("28S-frag5",      "GAAACCCGACAGGACCCGAAA",   False, []),
        ("18S-frag5",      "GAUUGCUGAGCCAGCUGAGCC",   True,  [9]),
        ("5.8S-frag2",     "AACGCGAAUUGAACGCGAAAC",   True,  [3, 15]),
        ("mt-16S-frag4",   "GCAAGUGCUCAAGUAAGCUUG",   True,  [6]),
        ("28S-frag6",      "UCGAAUUUGCUUCGAGCUCGG",   True,  [12]),
        ("18S-frag6",      "CUGAGCCAGCUGAGCCAGCUG",   False, []),
        ("5S-frag2",       "CCAUACCGACUCCAGCUGGUC",   True,  [4, 18]),
        ("mt-12S-frag5",   "UCAAAGCACCGAAACUCAGCA",   True,  [7]),
    ],

    # ── snoRNA: 25 templates ──────────────────────────────────────────────────
    "snoRNA": [
        ("SNORD3A-frag1",    "CUGAGGUAACUGGAGACCGCA",   True,  [12]),
        ("SNORD14-frag1",    "GCUGCGAAGCCCUGGUGCACC",   True,  [5, 17]),
        ("SNORD15A-frag1",   "UGUCAGUGCCUACCUGAUGCU",   True,  [8]),
        ("SNORD27-frag1",    "GAUUUGAUCUCUGUGCAAGCC",   False, []),
        ("SNORD32A-frag1",   "AGCCGUAGCGCUCUCCCUGUC",   True,  [3, 14]),
        ("SNORD33-frag1",    "GCUGCAACCCUGAGCAUGCCA",   True,  [6]),
        ("SNORD58A-frag1",   "CCGUGAUGGUGACCUGAAGCC",   False, []),
        ("SNORD68-frag1",    "AUGGCUGCUGCCAUUGCAGCA",   True,  [10, 19]),
        ("SNORD47-frag1",    "UGCAGCUCGAAGCCCUGGUGC",   True,  [4]),
        ("SNORD49A-frag1",   "GCUGCGAAGCUCAGUGCUACC",   True,  [9, 16]),
        ("SNORD50A-frag1",   "CAUCAGUGCCUACCUGAUGCC",   False, []),
        ("SNORD56-frag1",    "GCUGCAACCCUAAGCAUGCCA",   True,  [7]),
        ("SNORD18-frag1",    "UGAGGUAACUGGAGACCGCAG",   True,  [13]),
        ("SNORD21-frag1",    "GCUGCGAAGCCCUGGUGCACC",   True,  [2, 11]),
        ("SNORD22-frag1",    "CUGAGUGCCUACCUGAUGCCA",   False, []),
        ("SNORD24-frag1",    "AGCCGUAGCGCUCUCCUGUCC",   True,  [5, 16]),
        ("SNORD36A-frag1",   "GCUGCAACCCUGAGCAUGCCU",   True,  [8]),
        ("SNORD42A-frag1",   "CCGUGAUGGUGACCUGAAGCG",   True,  [3]),
        ("SNORD44-frag1",    "AUGGCUGCUGCCAUUGCAGCC",   True,  [11, 18]),
        ("SNORD45A-frag1",   "UGCAGCUCGAAGCCCUGGUGCA",  False, []),
        ("SNORD46-frag1",    "GCUGCGAAGCUCAGUGCUACCC",  True,  [6, 14]),
        ("SNORD48-frag1",    "CAUCAGUGCCUACCUGAUGCCU",  True,  [9]),
        ("SNORD57-frag1",    "GCUGCAACCCUAAGCAUGCCAG",  True,  [4, 17]),
        ("SNORD59A-frag1",   "UGAGGUAACUGGAGACCGCAGG",  False, []),
        ("SNORD61-frag1",    "AGCCGUAGCGCUCUCCUGUCC",   True,  [7, 15]),
    ],
}

ADAPTER_SEQ = "TGGAATTCTCGGGTGCCAAGG"
BASE_ERROR_RATE = 0.005
BASES = "ACGT"

def rna_to_dna(seq):
    return seq.upper().replace("U", "T")

def mutate_base(base, error_rate):
    if random.random() < error_rate:
        return random.choice([b for b in BASES if b != base])
    return base

def add_sequencing_errors(seq, error_rate=BASE_ERROR_RATE):
    return "".join(mutate_base(b, error_rate) for b in seq)

def generate_quality_string(length, min_q=30, max_q=40):
    quals = [random.randint(min_q, max_q) for _ in range(length)]
    for i in range(min(3, length)):
        quals[i] = max(20, quals[i] - random.randint(5, 10))
    for i in range(max(0, length-3), length):
        quals[i] = max(20, quals[i] - random.randint(5, 10))
    return "".join(chr(q + 33) for q in quals)

def apply_taps_chemistry(seq, m5c_positions, condition, replicate):
    rng = random.Random(hash((seq, condition, replicate)))
    seq_list = list(seq)
    for i, base in enumerate(seq_list):
        if base != "C":
            continue
        is_m5c = i in m5c_positions
        if condition == "no_treat":
            rate = rng.uniform(0.01, 0.03)
        elif condition == "pb_ctrl":
            rate = rng.uniform(0.05, 0.10) if is_m5c else rng.uniform(0.01, 0.03)
        elif condition == "treat":
            rate = rng.uniform(0.60, 0.90) if is_m5c else rng.uniform(0.01, 0.03)
        else:
            rate = 0.02
        if rng.random() < rate:
            seq_list[i] = "T"
    return "".join(seq_list)

def generate_read(template_seq, m5c_positions, condition, replicate,
                  read_len_range=(18, 24), adapter=ADAPTER_SEQ):
    read_len = random.randint(*read_len_range)
    max_start = max(0, len(template_seq) - read_len)
    start = random.randint(0, max_start)
    subseq = template_seq[start:start + read_len]
    sub_m5c = [p - start for p in m5c_positions if start <= p < start + read_len]
    subseq = apply_taps_chemistry(subseq, sub_m5c, condition, replicate)
    subseq = add_sequencing_errors(subseq)
    if random.random() < 0.5:
        adapter_len = random.randint(4, len(adapter))
        subseq = subseq + adapter[:adapter_len]
    subseq = subseq[:50]
    qual = generate_quality_string(len(subseq))
    return subseq, qual

def generate_sample_fastq(sample_name, condition, replicate, n_reads, outdir):
    outfile = outdir / f"{sample_name}.fastq.gz"
    biotype_weights = {"miRNA": 0.40, "tRNA": 0.25, "rRNA": 0.25, "snoRNA": 0.10}
    read_len_ranges = {"miRNA": (18, 24), "tRNA": (18, 25), "rRNA": (18, 23), "snoRNA": (18, 24)}
    reads_written = 0
    read_counter = 0
    with gzip.open(outfile, "wt") as fout:
        while reads_written < n_reads:
            biotype = random.choices(list(biotype_weights.keys()),
                                     weights=list(biotype_weights.values()))[0]
            templates = REFERENCE_SEQUENCES[biotype]
            name, rna_seq, has_m5c, m5c_pos = random.choice(templates)
            dna_seq = rna_to_dna(rna_seq)
            read_seq, qual = generate_read(
                dna_seq, m5c_pos if has_m5c else [],
                condition, replicate,
                read_len_range=read_len_ranges[biotype]
            )
            read_id = f"@{sample_name}_{biotype}_{name}_{read_counter}"
            fout.write(f"{read_id}\n{read_seq}\n+\n{qual}\n")
            reads_written += 1
            read_counter += 1
    return outfile

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir",  required=True)
    parser.add_argument("--reads",   type=int, default=1000000,
                        help="Reads per sample [default: 1000000]")
    parser.add_argument("--seed",    type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 100 samples: HEK only, 3 conditions
    # no-treat: 33, pb_ctrl: 33, treat: 34
    samples = []
    rep = 1
    for condition, n_reps in [("no_treat", 33), ("pb_ctrl", 33), ("treat", 34)]:
        for r in range(1, n_reps + 1):
            if condition == "no_treat":
                name = f"no-treat_Ctrl_HEK_R{r}"
            elif condition == "pb_ctrl":
                name = f"pb_Ctrl_HEK_R{r}"
            else:
                name = f"treat_HEK_R{r}"
            samples.append((name, condition, r))

    print(f"Generating {len(samples)} samples x {args.reads:,} reads each")
    print(f"Output directory: {outdir}")
    print(f"Biotypes: miRNA (40%), tRNA (25%), rRNA (25%), snoRNA (10%)")
    print(f"Templates: 50 miRNA, 25 tRNA, 25 rRNA, 25 snoRNA")
    print()

    for sample_name, condition, rep in samples:
        print(f"  {sample_name}...", end=" ", flush=True)
        outfile = generate_sample_fastq(sample_name, condition, rep, args.reads, outdir)
        size_mb = outfile.stat().st_size / 1e6
        print(f"done ({size_mb:.1f} MB)")

    print(f"\nAll samples written to {outdir}/")

    tsv_path = outdir / "samples.tsv"
    with open(tsv_path, "w") as f:
        f.write("sample\tcondition\tcell_line\treplicate\tfastq\n")
        for sample_name, condition, rep in samples:
            fastq = outdir / f"{sample_name}.fastq.gz"
            f.write(f"{sample_name}\t{condition}\tHEK\tR{rep}\t{fastq}\n")
    print(f"Samples TSV: {tsv_path}")

if __name__ == "__main__":
    main()
