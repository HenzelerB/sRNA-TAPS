#!/usr/bin/env python3
"""
simulate_taps_srna.py
=====================
Generate synthetic TAPS small RNA FASTQ files for sRNA-TAPS pipeline testing.

Simulates 9 samples (3 conditions x 3 replicates):
  - no_treat_HEK_R1/R2/R3  : background C→T only (2-5%)
  - pb_ctrl_HEK_R1/R2/R3   : background + slight PB signal (5-10% at m5C sites)
  - treat_HEK_R1/R2/R3     : genuine TAPS signal (40-80% C→T at m5C sites)

Chemistry: TAPS — modified cytosines appear as C→T in reads.
Adapter: TruSeq small RNA 3' adapter (TGGAATTCTCGGGTGCCAAGG)
Reference: Uses real hg38 sequences for miRNA/tRNA/rRNA/snoRNA loci.

Usage:
    python3 simulate_taps_srna.py --outdir /path/to/fastq_dir --reads 100000
"""

import argparse
import gzip
import os
import random
import sys
from pathlib import Path

# ── Seed sequences from real hg38 loci ───────────────────────────────────────
# Each entry: (name, sequence, biotype, has_m5C, m5C_positions)
# Sequences are representative mature/precursor sequences from miRBase/GtRNAdb
# m5C positions are 0-based within the sequence

REFERENCE_SEQUENCES = {

    # ── miRNA (mature, 18-23 nt) ─────────────────────────────────────────────
    "miRNA": [
        # hsa-miR-21-5p (chr17:59841266) — known m5C at pos 14
        ("hsa-miR-21-5p",   "UAGCUUAUCAGACUGAUGUUGA", True,  [14]),
        # hsa-miR-155-5p (chr21:26925472)
        ("hsa-miR-155-5p",  "UUAAUGCUAAUUGUGAUAGGGGU", True, [8, 18]),
        # hsa-miR-122-5p (chr18:56727408) — liver miRNA, m5C at C9
        ("hsa-miR-122-5p",  "UGGAGUGUGACAAUGGUGUUUG", True,  [9]),
        # hsa-miR-let-7a-5p (chr9:96938239)
        ("hsa-let-7a-5p",   "UGAGGUAGUAGGUUGUAUAGUU", False, []),
        # hsa-miR-16-5p (chr13:50048992)
        ("hsa-miR-16-5p",   "UAGCAGCACGUAAAUAUUGGCG", False, []),
        # hsa-miR-21-3p
        ("hsa-miR-21-3p",   "CAACACCAGUCGAUGGGCUGU", True,  [3, 15]),
        # hsa-miR-126-3p (chr9:136670613)
        ("hsa-miR-126-3p",  "UCGUACCGUGAGUAAUAAUGCG", True,  [2]),
        # hsa-miR-29a-3p
        ("hsa-miR-29a-3p",  "UAGCACCAUCUGAAAUCGGUUA", False, []),
        # hsa-miR-223-3p (chrX:65238686)
        ("hsa-miR-223-3p",  "UGUCAGUUUGUCAAAUACCCCA", True,  [19]),
        # hsa-miR-92a-3p
        ("hsa-miR-92a-3p",  "UAUUGCACUUGUCCCGGCCUGU", False, []),
    ],

    # ── tRNA fragments (mt-tRNA, 18-25 nt fragments) ─────────────────────────
    "tRNA": [
        # mt-tRNA-Leu(UUR) — chrM:3307, m5C at pos 34 (wobble, known m5C)
        ("mt-tRNA-Leu-frag1", "GGUUCGAUUCCCGGUCUCAGG", True,  [5, 12]),
        # mt-tRNA-Ile — chrM:4263
        ("mt-tRNA-Ile-frag1", "AGUGGUUUAAUCUUCUUAAGU", True,  [3, 17]),
        # mt-tRNA-Phe — chrM:577, m5C48 and m5C49 well characterised
        ("mt-tRNA-Phe-frag1", "GUCAUUCGUUGAAGAGUUGCA", True,  [8, 14]),
        # mt-tRNA-Val — chrM:1602
        ("mt-tRNA-Val-frag1", "CAUCAACUCCUAACACUUUCA", True,  [6]),
        # cytoplasmic tRNA-Gly fragment
        ("tRNA-Gly-frag1",    "GCGGUAGUAGCUCAGUCGGUA", False, []),
        # cytoplasmic tRNA-Ala fragment
        ("tRNA-Ala-frag1",    "GGGGGUAUAGCUCAGUGGUAG", True,  [4, 11]),
        # mt-tRNA-Pro — chrM:15956
        ("mt-tRNA-Pro-frag1", "CAAUCCAGUGCUUGAGUCACA", True,  [7]),
        # tRNA-Ser fragment
        ("tRNA-Ser-frag1",    "GCCCGGAUAGCUCAGUCGGUA", False, []),
    ],

    # ── rRNA fragments (18-25 nt fragments from 28S/16S/12S) ─────────────────
    "rRNA": [
        # 12S mt-rRNA — chrM:648, known m5C region
        ("mt-12S-frag1",  "CAUCACGAAACUCAGCACACU", True,  [4, 16]),
        ("mt-12S-frag2",  "UGGCUACACCUUGACAGCUAC", True,  [10]),
        # 16S mt-rRNA — chrM:1671
        ("mt-16S-frag1",  "GCUCGCCCUUGUGCAGAGAAU", True,  [2, 14]),
        ("mt-16S-frag2",  "AUGGCUGAGCCAGGCCUUUGA", True,  [8]),
        # 28S rRNA fragment (cytoplasmic)
        ("28S-frag1",     "GGCCCGAAACCCGACAGGACC", False, []),
        ("28S-frag2",     "CCGAGCUCGAAUUUGCUUCGA", True,  [5]),
        # 18S rRNA fragment
        ("18S-frag1",     "ACGGUAGAGCUACCGAUUGCU", False, []),
        ("18S-frag2",     "CGAAACUCGCCCAGCAAACGC", True,  [3, 18]),
        # 5.8S rRNA
        ("5.8S-frag1",    "GCGAAACGCGAAUUGAACGCG", True,  [9]),
    ],

    # ── snoRNA fragments (18-25 nt) ───────────────────────────────────────────
    "snoRNA": [
        # SNORD3A (RNU3-1) — chr17
        ("SNORD3A-frag1",  "CUGAGGUAACUGGAGACCGCA", True,  [12]),
        # SNORD14 (SNORD14A) — chr1
        ("SNORD14-frag1",  "GCUGCGAAGCCCUGGUGCACC", True,  [5, 17]),
        # SNORD15A — chr11
        ("SNORD15A-frag1", "UGUCAGUGCCUACCUGAUGCU", True,  [8]),
        # SNORD27 — chr1
        ("SNORD27-frag1",  "GAUUUGAUCUCUGUGCAAGCC", False, []),
        # SNORD32A
        ("SNORD32A-frag1", "AGCCGUAGCGCUCUCCCUGUC", True,  [3, 14]),
        # SNORD33
        ("SNORD33-frag1",  "GCUGCAACCCUGAGCAUGCCA", True,  [6]),
        # SNORD58A
        ("SNORD58A-frag1", "CCGUGAUGGUGACCUGAAGCC", False, []),
        # SNORD68
        ("SNORD68-frag1",  "AUGGCUGCUGCCAUUGCAGCA", True,  [10, 19]),
    ],
}

# ── TruSeq small RNA adapter ──────────────────────────────────────────────────
ADAPTER_SEQ = "TGGAATTCTCGGGTGCCAAGG"

# ── Sequencing error profile ──────────────────────────────────────────────────
BASE_ERROR_RATE = 0.005   # 0.5% per base
BASES = "ACGT"

# RNA→DNA conversion (U→T)
def rna_to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")

def complement(base: str) -> str:
    return {"A":"T","T":"A","G":"C","C":"G","N":"N"}.get(base.upper(), "N")

def reverse_complement(seq: str) -> str:
    return "".join(complement(b) for b in reversed(seq))

def mutate_base(base: str, error_rate: float) -> str:
    if random.random() < error_rate:
        return random.choice([b for b in BASES if b != base])
    return base

def add_sequencing_errors(seq: str, error_rate: float = BASE_ERROR_RATE) -> str:
    return "".join(mutate_base(b, error_rate) for b in seq)

def generate_quality_string(length: int, min_q: int = 30, max_q: int = 40) -> str:
    """Generate Phred+33 quality string."""
    quals = [random.randint(min_q, max_q) for _ in range(length)]
    # Lower quality at ends (realistic)
    for i in range(min(3, length)):
        quals[i] = max(20, quals[i] - random.randint(5, 10))
    for i in range(max(0, length-3), length):
        quals[i] = max(20, quals[i] - random.randint(5, 10))
    return "".join(chr(q + 33) for q in quals)

def apply_taps_chemistry(seq: str, m5c_positions: list,
                         condition: str, replicate: int,
                         background_rate: float = 0.03) -> str:
    """
    Apply TAPS chemistry: modified C → T.
    
    no_treat : background C→T only (2-5%)
    pb_ctrl  : background + slight signal at m5C positions (5-15%)
    treat    : genuine m5C signal (40-80% C→T) + background elsewhere
    """
    rng = random.Random(hash((seq, condition, replicate)))
    seq_list = list(seq)
    
    for i, base in enumerate(seq_list):
        if base != "C":
            continue
        
        is_m5c = i in m5c_positions
        
        if condition == "no_treat":
            rate = rng.uniform(0.02, 0.05)
        elif condition == "pb_ctrl":
            rate = rng.uniform(0.05, 0.15) if is_m5c else rng.uniform(0.02, 0.05)
        elif condition == "treat":
            rate = rng.uniform(0.40, 0.80) if is_m5c else rng.uniform(0.02, 0.05)
        else:
            rate = 0.03
        
        if rng.random() < rate:
            seq_list[i] = "T"
    
    return "".join(seq_list)

def generate_read(template_seq: str, m5c_positions: list,
                  condition: str, replicate: int,
                  read_len_range: tuple = (18, 24),
                  adapter: str = ADAPTER_SEQ) -> tuple:
    """
    Generate one synthetic read from a template sequence.
    Returns (read_sequence, quality_string).
    """
    # Random read length within biotype range
    read_len = random.randint(*read_len_range)
    
    # Random start within template (allow slight trimming)
    max_start = max(0, len(template_seq) - read_len)
    start = random.randint(0, max_start)
    
    # Extract subsequence
    subseq = template_seq[start:start + read_len]
    
    # Adjust m5C positions to subsequence coordinates
    sub_m5c = [p - start for p in m5c_positions
               if start <= p < start + read_len]
    
    # Apply TAPS chemistry
    subseq = apply_taps_chemistry(subseq, sub_m5c, condition, replicate)
    
    # Add sequencing errors
    subseq = add_sequencing_errors(subseq)
    
    # Add partial adapter (50% of reads have adapter, realistic for short fragments)
    if random.random() < 0.5:
        adapter_len = random.randint(4, len(adapter))
        subseq = subseq + adapter[:adapter_len]
    
    # Truncate to max read length (typical Illumina 50bp)
    subseq = subseq[:50]
    
    qual = generate_quality_string(len(subseq))
    return subseq, qual

def generate_sample_fastq(sample_name: str, condition: str, replicate: int,
                           n_reads: int, outdir: Path) -> Path:
    """Generate one FASTQ.gz file for a sample."""
    
    outfile = outdir / f"{sample_name}.fastq.gz"
    
    # Biotype proportions (realistic for small RNA-seq)
    biotype_weights = {
        "miRNA":  0.40,
        "tRNA":   0.25,
        "rRNA":   0.25,
        "snoRNA": 0.10,
    }
    
    # Read length ranges per biotype
    read_len_ranges = {
        "miRNA":  (18, 24),
        "tRNA":   (18, 25),
        "rRNA":   (18, 23),
        "snoRNA": (18, 24),
    }
    
    reads_written = 0
    read_counter = 0
    
    with gzip.open(outfile, "wt") as fout:
        while reads_written < n_reads:
            # Pick a biotype
            biotype = random.choices(
                list(biotype_weights.keys()),
                weights=list(biotype_weights.values())
            )[0]
            
            # Pick a template from that biotype
            templates = REFERENCE_SEQUENCES[biotype]
            name, rna_seq, has_m5c, m5c_pos = random.choice(templates)
            
            # Convert RNA to DNA
            dna_seq = rna_to_dna(rna_seq)
            
            # Generate read
            read_seq, qual = generate_read(
                dna_seq, m5c_pos if has_m5c else [],
                condition, replicate,
                read_len_range=read_len_ranges[biotype]
            )
            
            # Write FASTQ entry
            read_id = f"@{sample_name}_{biotype}_{name}_{read_counter}"
            fout.write(f"{read_id}\n{read_seq}\n+\n{qual}\n")
            
            reads_written += 1
            read_counter += 1
    
    return outfile

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir",  required=True,
                        help="Output directory for FASTQ files")
    parser.add_argument("--reads",   type=int, default=100000,
                        help="Number of reads per sample [default: 100000]")
    parser.add_argument("--seed",    type=int, default=42,
                        help="Random seed for reproducibility [default: 42]")
    args = parser.parse_args()
    
    random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Define 9 samples: 3 conditions x 3 replicates, single cell line HEK
    samples = []
    for condition in ["no_treat", "pb_ctrl", "treat"]:
        for rep in [1, 2, 3]:
            # Match your actual sample naming convention
            if condition == "no_treat":
                name = f"no-treat_Ctrl_HEK_R{rep}"
            elif condition == "pb_ctrl":
                name = f"pb_Ctrl_HEK_R{rep}"
            else:
                name = f"treat_HEK_R{rep}"
            samples.append((name, condition, rep))
    
    print(f"Generating {len(samples)} samples x {args.reads:,} reads each")
    print(f"Output directory: {outdir}")
    print(f"Biotypes: miRNA, tRNA, rRNA, snoRNA")
    print()
    
    for sample_name, condition, rep in samples:
        print(f"  Generating {sample_name} ({condition}, R{rep})...", end=" ", flush=True)
        outfile = generate_sample_fastq(
            sample_name, condition, rep,
            args.reads, outdir
        )
        size_mb = outfile.stat().st_size / 1e6
        print(f"done ({size_mb:.1f} MB)")
    
    print(f"\nAll samples written to {outdir}/")
    print("\nSample sheet (samples.tsv):")
    print("sample\tcondition\tcell_line\treplicate\tfastq")
    for sample_name, condition, rep in samples:
        fastq = outdir / f"{sample_name}.fastq.gz"
        print(f"{sample_name}\t{condition}\tHEK\tR{rep}\t{fastq}")
    
    # Write samples.tsv
    tsv_path = outdir / "samples.tsv"
    with open(tsv_path, "w") as f:
        f.write("sample\tcondition\tcell_line\treplicate\tfastq\n")
        for sample_name, condition, rep in samples:
            fastq = outdir / f"{sample_name}.fastq.gz"
            f.write(f"{sample_name}\t{condition}\tHEK\tR{rep}\t{fastq}\n")
    print(f"\nSamples TSV: {tsv_path}")


if __name__ == "__main__":
    main()
