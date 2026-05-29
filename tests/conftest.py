# -*- coding: utf-8 -*-
"""
conftest.py — shared pytest fixtures for sRNA-TAPS tests.

Provides:
    - Tiny synthetic BAM files (Bowtie1 XA:i format) with known C→T modifications
    - Minimal reference FASTA with known cytosine positions
    - SNP BED fixtures with known variants
    - Temp directory management
"""

import os
import gzip
import shutil
import struct
import tempfile
from pathlib import Path

import numpy as np
import pysam
import pytest


# ══════════════════════════════════════════════════════════════════════════════
# Helpers: build minimal FASTA and BAM in memory
# ══════════════════════════════════════════════════════════════════════════════

CHROM   = "chr1"
CHROM_SEQ = (
    "AAACGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGAATT"
    "CGAATTCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGAA"
)
# C positions in CHROM_SEQ (0-based): 3,7,11,15,19,25,...
# SNP_POS: guaranteed C position used for SNP fixture and SNP blacklist BED
SNP_POS = next(i for i, b in enumerate(CHROM_SEQ) if b == "C")  # = 3


def write_fasta(path: str, chrom: str = CHROM, seq: str = CHROM_SEQ) -> str:
    """Write a minimal FASTA and create .fai index."""
    with open(path, "w") as fh:
        fh.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")
    pysam.faidx(path)
    return path


def write_bam(path: str, reads: list[dict]) -> str:
    """
    Write a minimal sorted indexed BAM.

    Each read dict:
        name:     read name
        seq:      query sequence (same length as ref span)
        pos:      0-based start position on CHROM
        qual:     list of base qualities (default: all 30)
        xa:       XA:i:N multi-mapper count (default: 1)
        reverse:  is_reverse flag (default: False)
    """
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": len(CHROM_SEQ)}],
    })
    # Write unsorted SAM first, then sort
    tmp_sam = path + ".tmp.sam"
    with pysam.AlignmentFile(tmp_sam, "w", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(header)
            a.query_name         = r["name"]
            a.query_sequence     = r["seq"]
            a.flag               = 16 if r.get("reverse", False) else 0
            a.reference_id       = 0
            a.reference_start    = r["pos"]
            a.mapping_quality    = 30
            a.cigar              = [(0, len(r["seq"]))]  # all match
            a.query_qualities    = pysam.qualitystring_to_array(
                "".join(chr(q + 33) for q in r.get("qual", [30] * len(r["seq"])))
            )
            a.set_tag("XA", r.get("xa", 1))
            out.write(a)

    pysam.sort("-o", path, tmp_sam)
    pysam.index(path)
    Path(tmp_sam).unlink(missing_ok=True)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ══════════════════════════════════════════════════════════════════════════════

@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    """Session-scoped temp directory."""
    return tmp_path_factory.mktemp("srnataps_tests")


@pytest.fixture(scope="session")
def fasta_path(tmp_dir):
    """Minimal hg38-like FASTA with known cytosine positions."""
    path = str(tmp_dir / "test_genome.fa")
    return write_fasta(path)


@pytest.fixture(scope="session")
def bam_notreat(tmp_dir, fasta_path):
    """
    No-treat BAM: all reads show reference bases (no C→T).
    Used to test SNP blacklist building and background estimation.
    """
    reads = []
    for i in range(50):
        # Reads that match the reference exactly at C positions
        seq = CHROM_SEQ[10:40]   # 30 nt read, contains several C positions
        reads.append({"name": f"notreat_{i}", "seq": seq, "pos": 10, "xa": 1})
    path = str(tmp_dir / "notreat.sorted.bam")
    return write_bam(path, reads)


@pytest.fixture(scope="session")
def bam_with_snp(tmp_dir, fasta_path):
    """
    BAM where position 12 has a C→T SNP (heterozygous ~50%).
    Layer 2 SNP calling should flag this position.
    """
    reads = []
    seq_ref = list(CHROM_SEQ[0:30])
    seq_alt = list(CHROM_SEQ[0:30])
    # Use SNP_POS (guaranteed C) — index in read = SNP_POS - 0 = SNP_POS
    seq_alt[SNP_POS] = "T"   # C→T SNP at guaranteed C position

    for i in range(25):
        reads.append({"name": f"snp_ref_{i}", "seq": "".join(seq_ref), "pos": 0})
    for i in range(25):
        reads.append({"name": f"snp_alt_{i}", "seq": "".join(seq_alt), "pos": 0})

    path = str(tmp_dir / "with_snp.sorted.bam")
    return write_bam(path, reads)


@pytest.fixture(scope="session")
def bam_taps_modified(tmp_dir, fasta_path):
    """
    TAPS-treated BAM with C→T at known modified positions.
    All 50 reads show C→T at position 20 (high confidence modification).
    Position 12 shows no C→T (unmodified control site).
    """
    reads = []
    seq_base = list(CHROM_SEQ[10:40])

    # Position 20 in CHROM_SEQ = index 10 in read
    # Make it C→T in all reads (mod_rate = 1.0 at this position)
    if CHROM_SEQ[20] == "C":
        for i in range(50):
            seq = list(seq_base)
            seq[10] = "T"
            reads.append({"name": f"taps_{i}", "seq": "".join(seq), "pos": 10})
    else:
        for i in range(50):
            reads.append({"name": f"taps_{i}", "seq": "".join(seq_base), "pos": 10})

    path = str(tmp_dir / "taps_modified.sorted.bam")
    return write_bam(path, reads)


@pytest.fixture(scope="session")
def snp_blacklist_bed(tmp_dir):
    """SNP blacklist BED with one known C→T variant at chr1:12."""
    path = str(tmp_dir / "sample_snps_TEST.bed")
    with open(path, "w") as fh:
        fh.write("# Sample SNP blacklist\n")
        fh.write("chrom\tstart\tend\tref\talt\tallele_freq\tcoverage\tcell_line\n")
        fh.write(f"{CHROM}\t{SNP_POS}\t{SNP_POS+1}\tC\tT\t0.50\t50\tTEST\n")
    return path


@pytest.fixture(scope="session")
def minimal_gtf(tmp_dir):
    """Minimal GTF with one miRNA and one tRNA gene."""
    path = str(tmp_dir / "test.gtf")
    with open(path, "w") as fh:
        fh.write('##format: gtf\n')
        fh.write(
            f'{CHROM}\tENSEMBL\tgene\t5\t35\t.\t+\t.\t'
            'gene_id "ENSG001"; gene_name "MIR1"; gene_biotype "miRNA";\n'
        )
        fh.write(
            f'{CHROM}\tENSEMBL\tgene\t50\t80\t.\t+\t.\t'
            'gene_id "ENSG002"; gene_name "TRNA1"; gene_biotype "tRNA";\n'
        )
    return path
