# -*- coding: utf-8 -*-
"""
tests/test_caller.py — Tests for srnataps.caller

Tests:
    - process_chromosome: correct C→T counting with multi-mapper weighting
    - SNP sites excluded from pileup (correct ordering — before counting)
    - binomial_test_bh: correct p-value and FDR correction
    - Full calling pipeline: PASS sites correctly identified
    - SNP_flag column present and correct in output TSV
"""

import os
import sys
import pandas as pd
import pytest
from pathlib import Path

from srnataps.caller import (
    PASS, SNP_KNOWN, SNP_SAMPLE, SNP_HET, SNP_MULTI,
    process_chromosome,
    binomial_test_bh,
    load_dbsnp,
    load_sample_snp_bed,
)
from conftest import CHROM, CHROM_SEQ, SNP_POS


# ══════════════════════════════════════════════════════════════════════════════
# process_chromosome
# ══════════════════════════════════════════════════════════════════════════════

class TestProcessChromosome:

    def test_counts_ct_modification(self, bam_taps_modified, fasta_path):
        """
        BAM has C→T at position 20. process_chromosome should count T at that position.
        """
        job = (
            bam_taps_modified, fasta_path, CHROM,
            20, 10, "ALL",
            frozenset(), frozenset(), frozenset(),  # no SNP filters
        )
        rows = process_chromosome(job)

        # Find position 20
        pos_snp = [r for r in rows if r["start"] == SNP_POS]

        if CHROM_SEQ[20] == "C":
            assert len(pos20) == 1, "Position 20 should be called (coverage >= 10)"
            assert pos20[0]["mod_count"] > 0, "Should have mod_count > 0"
            assert pos20[0]["mod_rate"] > 0.8, "mod_rate should be high (all reads modified)"

    def test_snp_site_excluded_before_counting(self, bam_taps_modified, fasta_path):
        """
        Critical test: SNP at position 20 should result in ZERO counts,
        not a site with mod_count > 0 that gets filtered later.
        SNP filtering must happen BEFORE accumulating counts.
        """
        # Flag position 20 as a known SNP
        snp_key = f"{CHROM}:{SNP_POS}"
        dbsnp = frozenset([snp_key])

        job = (
            bam_taps_modified, fasta_path, CHROM,
            20, 10, "ALL",
            dbsnp, frozenset(), frozenset(),
        )
        rows = process_chromosome(job)

        # Position 20 should NOT appear in the output at all
        pos_snp = [r for r in rows if r["start"] == SNP_POS]
        assert len(pos_snp) == 0, (
            "SNP-flagged site should be excluded before counting — "
            "not counted and then filtered. If this fails, SNP filtering "
            "is happening after counting (incorrect order)."
        )

    def test_multimapper_weighting(self, tmp_dir, fasta_path):
        """
        Reads with XA:i:2 should contribute weight 0.5 each.
        Two such reads at a C position → mod_count = 1.0, not 2.0.
        """
        from conftest import write_bam, CHROM_SEQ

        # Find a C position to modify
        c_pos = next(i for i, b in enumerate(CHROM_SEQ[15:35], start=15) if b == "C")

        reads = []
        for i in range(20):
            seq = list(CHROM_SEQ[10:40])
            seq[c_pos - 10] = "T"
            reads.append({
                "name": f"multi_{i}",
                "seq":  "".join(seq),
                "pos":  10,
                "xa":   2,   # maps to 2 loci → weight = 0.5
            })

        bam = write_bam(str(tmp_dir / "multi_test.sorted.bam"), reads)

        job = (bam, fasta_path, CHROM, 20, 5, "ALL",
               frozenset(), frozenset(), frozenset())
        rows = process_chromosome(job)

        pos = [r for r in rows if r["start"] == c_pos]
        if pos:
            # 20 reads × weight 0.5 = 10.0, not 20.0
            assert abs(pos[0]["mod_count"] - 10.0) < 0.1, (
                f"Expected mod_count ~10.0 with XA=2, got {pos[0]['mod_count']}"
            )

    def test_min_coverage_filter(self, bam_notreat, fasta_path):
        """Sites with coverage < min_cov should not appear in output."""
        job = (
            bam_notreat, fasta_path, CHROM,
            20, 1000,  # extremely high min_cov — nothing should pass
            "ALL",
            frozenset(), frozenset(), frozenset(),
        )
        rows = process_chromosome(job)
        assert len(rows) == 0, "No sites should pass with min_cov=1000"

    def test_output_columns_present(self, bam_taps_modified, fasta_path):
        """Output rows should have all required columns."""
        job = (
            bam_taps_modified, fasta_path, CHROM,
            20, 5, "ALL",
            frozenset(), frozenset(), frozenset(),
        )
        rows = process_chromosome(job)
        if rows:
            required = {"chrom", "start", "end", "context",
                        "mod_count", "unmod_count", "coverage", "mod_rate", "snp_flag"}
            assert required.issubset(set(rows[0].keys()))


# ══════════════════════════════════════════════════════════════════════════════
# binomial_test_bh
# ══════════════════════════════════════════════════════════════════════════════

class TestBinomialTestBH:

    def test_high_mod_site_is_significant(self):
        """A site with mod_rate = 1.0 at high coverage should be highly significant."""
        rows = [{"mod_count": 50.0, "unmod_count": 0.0, "coverage": 50.0,
                 "mod_rate": 1.0, "snp_flag": PASS}]
        result = binomial_test_bh(rows, background_rate=0.005)
        assert result[0]["pvalue"] < 1e-10
        assert result[0]["padj"]   < 0.05

    def test_low_mod_site_is_not_significant(self):
        """A site at background rate should not be significant."""
        rows = [{"mod_count": 1.0, "unmod_count": 199.0, "coverage": 200.0,
                 "mod_rate": 0.005, "snp_flag": PASS}]
        result = binomial_test_bh(rows, background_rate=0.005)
        assert result[0]["pvalue"] > 0.05

    def test_bh_correction_applied(self):
        """BH correction should increase p-values when many tests are run."""
        rows = [
            {"mod_count": 50.0, "unmod_count": 0.0,   "coverage": 50.0,  "mod_rate": 1.0,   "snp_flag": PASS},
            {"mod_count": 1.0,  "unmod_count": 199.0, "coverage": 200.0, "mod_rate": 0.005, "snp_flag": PASS},
        ]
        result = binomial_test_bh(rows, background_rate=0.005)
        # padj should be >= pvalue (BH never decreases p-values)
        for r in result:
            assert r["padj"] >= r["pvalue"] - 1e-10

    def test_empty_input_returns_empty(self):
        assert binomial_test_bh([], background_rate=0.005) == []

    def test_pvalue_and_padj_keys_added(self):
        rows = [{"mod_count": 10.0, "unmod_count": 0.0, "coverage": 10.0,
                 "mod_rate": 1.0, "snp_flag": PASS}]
        result = binomial_test_bh(rows, background_rate=0.005)
        assert "pvalue" in result[0]
        assert "padj"   in result[0]


# ══════════════════════════════════════════════════════════════════════════════
# Full calling pipeline (integration)
# ══════════════════════════════════════════════════════════════════════════════

class TestFullCalling:

    def test_output_tsv_created(self, bam_taps_modified, fasta_path,
                                snp_blacklist_bed, tmp_dir):
        """Full pipeline should produce a TSV with correct columns."""
        out_tsv = str(tmp_dir / "test_taps_output.tsv")

        sys.argv = [
            "07_taps_calling.py",
            "--bam",             bam_taps_modified,
            "--fasta",           fasta_path,
            "--out",             out_tsv,
            "--min-cov",         "5",
            "--min-qual",        "20",
            "--background-rate", "0.005",
            "--sample-snp-bed",  snp_blacklist_bed,
            "--het-threshold",   "0.40",
            "--cell-line",       "TEST",
            "--threads",         "1",
        ]
        from srnataps.caller import main as caller_main
        caller_main()

        assert Path(out_tsv).exists(), "Output TSV should be created"
        df = pd.read_csv(out_tsv, sep="\t")

        required_cols = {"chrom", "start", "end", "context",
                         "mod_count", "unmod_count", "coverage",
                         "mod_rate", "pvalue", "padj", "snp_flag"}
        assert required_cols.issubset(set(df.columns))

    def test_output_contains_only_pass_sites(self, bam_taps_modified, fasta_path,
                                              snp_blacklist_bed, tmp_dir):
        """All sites in output TSV should have snp_flag == PASS."""
        out_tsv = str(tmp_dir / "test_pass_only.tsv")

        sys.argv = [
            "07_taps_calling.py",
            "--bam",             bam_taps_modified,
            "--fasta",           fasta_path,
            "--out",             out_tsv,
            "--min-cov",         "5",
            "--min-qual",        "20",
            "--background-rate", "0.005",
            "--sample-snp-bed",  snp_blacklist_bed,
            "--het-threshold",   "0.40",
            "--cell-line",       "TEST",
            "--threads",         "1",
        ]
        from srnataps.caller import main as caller_main
        caller_main()

        df = pd.read_csv(out_tsv, sep="\t")
        if len(df) > 0:
            assert (df["snp_flag"] == PASS).all(), \
                "All output sites should have snp_flag == PASS"

    def test_snp_position_excluded(self, bam_taps_modified, fasta_path,
                                    snp_blacklist_bed, tmp_dir):
        """Position 12 (in SNP blacklist) should not appear in output."""
        out_tsv = str(tmp_dir / "test_no_snp_pos.tsv")

        sys.argv = [
            "07_taps_calling.py",
            "--bam",             bam_taps_modified,
            "--fasta",           fasta_path,
            "--out",             out_tsv,
            "--min-cov",         "5",
            "--min-qual",        "20",
            "--background-rate", "0.005",
            "--sample-snp-bed",  snp_blacklist_bed,
            "--het-threshold",   "0.40",
            "--cell-line",       "TEST",
            "--threads",         "1",
        ]
        from srnataps.caller import main as caller_main
        caller_main()

        df = pd.read_csv(out_tsv, sep="\t")
        snp_sites = df[(df["chrom"] == CHROM) & (df["start"] == SNP_POS)]
        assert len(snp_sites) == 0, \
            "SNP position 12 should not appear in PASS output"
