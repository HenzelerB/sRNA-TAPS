# -*- coding: utf-8 -*-
"""
tests/test_snp.py — Tests for srnataps.snp

Tests:
    - build_sample_snp_blacklist: correct C→T detection at known SNP position
    - load_sample_snp_bed: correct loading of positions and het flagging
    - assign_snp_flag: correct flag assignment for all three layers and combinations
"""

import pytest
from pathlib import Path

from srnataps.caller import (
    PASS, SNP_KNOWN, SNP_SAMPLE, SNP_HET, SNP_MULTI,
    assign_snp_flag,
    load_sample_snp_bed,
    load_dbsnp,
)
from srnataps.snp import main as snp_main
from conftest import CHROM, SNP_POS


# ══════════════════════════════════════════════════════════════════════════════
# assign_snp_flag
# ══════════════════════════════════════════════════════════════════════════════

class TestAssignSnpFlag:
    """Unit tests for the SNP flag assignment logic."""

    def test_pass_when_no_flags(self):
        key = "chr1:100"
        assert assign_snp_flag(key, frozenset(), frozenset(), frozenset()) == PASS

    def test_single_dbsnp(self):
        key = "chr1:100"
        assert assign_snp_flag(key, frozenset([key]), frozenset(), frozenset()) == SNP_KNOWN

    def test_single_sample(self):
        key = "chr1:100"
        assert assign_snp_flag(key, frozenset(), frozenset([key]), frozenset()) == SNP_SAMPLE

    def test_single_het(self):
        key = "chr1:100"
        assert assign_snp_flag(key, frozenset(), frozenset(), frozenset([key])) == SNP_HET

    def test_multiple_flags_returns_multi(self):
        key = "chr1:100"
        assert assign_snp_flag(
            key,
            frozenset([key]),
            frozenset([key]),
            frozenset(),
        ) == SNP_MULTI

    def test_all_three_flags_returns_multi(self):
        key = "chr1:100"
        assert assign_snp_flag(
            key,
            frozenset([key]),
            frozenset([key]),
            frozenset([key]),
        ) == SNP_MULTI

    def test_different_position_is_pass(self):
        key     = "chr1:100"
        snp_key = "chr1:200"
        assert assign_snp_flag(key, frozenset([snp_key]), frozenset(), frozenset()) == PASS

    def test_chrom_specificity(self):
        """Flag at chr1:100 should not affect chr2:100."""
        key1 = "chr1:100"
        key2 = "chr2:100"
        assert assign_snp_flag(key2, frozenset([key1]), frozenset(), frozenset()) == PASS


# ══════════════════════════════════════════════════════════════════════════════
# load_sample_snp_bed
# ══════════════════════════════════════════════════════════════════════════════

class TestLoadSampleSnpBed:

    def test_loads_known_position(self, snp_blacklist_bed):
        sample_snps, het_positions = load_sample_snp_bed(snp_blacklist_bed, het_threshold=0.40)
        assert f"{CHROM}:{SNP_POS}" in sample_snps

    def test_het_threshold_applied(self, snp_blacklist_bed):
        """Position with AF=0.50 should be in het_positions at threshold 0.40."""
        _, het_positions = load_sample_snp_bed(snp_blacklist_bed, het_threshold=0.40)
        assert f"{CHROM}:{SNP_POS}" in het_positions

    def test_het_threshold_not_applied_when_below(self, snp_blacklist_bed):
        """Position with AF=0.50 should NOT be in het_positions at threshold 0.60."""
        _, het_positions = load_sample_snp_bed(snp_blacklist_bed, het_threshold=0.60)
        assert f"{CHROM}:{SNP_POS}" not in het_positions

    def test_missing_file_returns_empty(self, tmp_dir):
        sample_snps, het = load_sample_snp_bed(
            str(tmp_dir / "nonexistent.bed"), het_threshold=0.40
        )
        assert len(sample_snps) == 0
        assert len(het) == 0

    def test_returns_frozensets(self, snp_blacklist_bed):
        sample_snps, het = load_sample_snp_bed(snp_blacklist_bed, het_threshold=0.40)
        assert isinstance(sample_snps, frozenset)
        assert isinstance(het, frozenset)


# ══════════════════════════════════════════════════════════════════════════════
# build_sample_snp_blacklist (integration)
# ══════════════════════════════════════════════════════════════════════════════

class TestBuildSnpBlacklist:

    def test_detects_ct_snp(self, bam_with_snp, fasta_path, tmp_dir):
        """
        BAM has C→T at position 12 at ~50% AF.
        build_sample_snp_blacklist should detect it.
        """
        from srnataps.snp import main as run_snp
        out_bed = str(tmp_dir / "snps_test.bed")

        import sys
        sys.argv = [
            "06_build_snp_blacklist.py",
            "--bam",       bam_with_snp,
            "--fasta",     fasta_path,
            "--out",       out_bed,
            "--min-af",    "0.20",
            "--min-cov",   "10",
            "--cell-line", "TEST",
        ]
        run_snp()

        assert Path(out_bed).exists()
        with open(out_bed) as fh:
            lines = [l for l in fh if not l.startswith("#") and "chrom" not in l]

        # Should have at least one SNP entry
        assert len(lines) >= 1

        # The SNP at position 12 should be detected
        positions = [l.split("\t")[1] for l in lines]
        assert str(SNP_POS) in positions

    def test_notreat_bam_has_no_snps(self, bam_notreat, fasta_path, tmp_dir):
        """
        No-treat BAM with no C→T transitions should produce empty blacklist.
        """
        from srnataps.snp import main as run_snp
        out_bed = str(tmp_dir / "snps_empty.bed")

        import sys
        sys.argv = [
            "06_build_snp_blacklist.py",
            "--bam",       bam_notreat,
            "--fasta",     fasta_path,
            "--out",       out_bed,
            "--min-af",    "0.20",
            "--min-cov",   "10",
            "--cell-line", "TEST",
        ]
        run_snp()

        with open(out_bed) as fh:
            lines = [l for l in fh if not l.startswith("#") and "chrom" not in l]
        assert len(lines) == 0
