# -*- coding: utf-8 -*-
"""
tests/test_annotate.py — Tests for srnataps.annotate

Tests:
    - parse_gtf: correct loading and interval building
    - find_biotype: correct priority assignment
    - build_chrom_aliases: chr prefix handling
    - Full biotype split: reads assigned to correct biotypes
"""

import pytest
from srnataps.annotate import (
    parse_gtf,
    find_biotype,
    build_chrom_aliases,
    BIOTYPE_PRIORITY,
    BIOTYPE_MAP,
)
from conftest import CHROM


class TestParseGtf:

    def test_loads_gene_features(self, minimal_gtf):
        intervals, starts = parse_gtf(minimal_gtf)
        assert len(intervals) > 0

    def test_miRNA_feature_loaded(self, minimal_gtf):
        intervals, starts = parse_gtf(minimal_gtf)
        assert CHROM in intervals
        biotypes = [feat[2] for feat in intervals[CHROM]]
        assert "miRNA" in biotypes

    def test_trna_feature_loaded(self, minimal_gtf):
        intervals, starts = parse_gtf(minimal_gtf)
        biotypes = [feat[2] for feat in intervals[CHROM]]
        assert "tRNA" in biotypes

    def test_starts_sorted(self, minimal_gtf):
        intervals, starts = parse_gtf(minimal_gtf)
        for chrom, start_list in starts.items():
            assert start_list == sorted(start_list), \
                f"Starts not sorted for {chrom}"


class TestFindBiotype:

    def test_miRNA_position_returns_mirna(self, minimal_gtf):
        """Position 15 is inside the miRNA gene (5-35)."""
        intervals, starts = parse_gtf(minimal_gtf)
        aliases = {CHROM: [CHROM]}
        result = find_biotype(CHROM, 15, intervals, starts, aliases)
        assert result == "miRNA"

    def test_trna_position_returns_trna(self, minimal_gtf):
        """Position 60 is inside the tRNA gene (50-80)."""
        intervals, starts = parse_gtf(minimal_gtf)
        aliases = {CHROM: [CHROM]}
        result = find_biotype(CHROM, 60, intervals, starts, aliases)
        assert result == "tRNA"

    def test_intergenic_returns_other(self, minimal_gtf):
        """Position 200 is outside any annotated gene."""
        intervals, starts = parse_gtf(minimal_gtf)
        aliases = {CHROM: [CHROM]}
        result = find_biotype(CHROM, 200, intervals, starts, aliases)
        assert result == "other"

    def test_priority_miRNA_over_tRNA(self, tmp_dir):
        """When miRNA and tRNA overlap, miRNA should win (priority 1 vs 2)."""
        import tempfile
        gtf_path = str(tmp_dir / "overlap_test.gtf")
        with open(gtf_path, "w") as fh:
            fh.write(
                f'{CHROM}\ttest\tgene\t10\t50\t.\t+\t.\t'
                'gene_id "G1"; gene_biotype "miRNA";\n'
            )
            fh.write(
                f'{CHROM}\ttest\tgene\t10\t50\t.\t+\t.\t'
                'gene_id "G2"; gene_biotype "tRNA";\n'
            )
        intervals, starts = parse_gtf(gtf_path)
        aliases = {CHROM: [CHROM]}
        result = find_biotype(CHROM, 25, intervals, starts, aliases)
        assert result == "miRNA", "miRNA should have priority over tRNA"

    def test_unknown_chrom_returns_other(self, minimal_gtf):
        intervals, starts = parse_gtf(minimal_gtf)
        aliases = {}
        result = find_biotype("chrUNKNOWN", 15, intervals, starts, aliases)
        assert result == "other"


class TestBuildChromAliases:

    def test_chr_prefix_handled(self):
        """BAM chrom 'chr1' should map to GTF chrom '1'."""
        aliases = build_chrom_aliases(["chr1"], ["1", "2", "3"])
        assert "1" in aliases.get("chr1", [])

    def test_no_prefix_to_chr(self):
        """BAM chrom '1' should map to GTF chrom 'chr1'."""
        aliases = build_chrom_aliases(["1"], ["chr1", "chr2"])
        assert "chr1" in aliases.get("1", [])

    def test_chrM_to_MT(self):
        """BAM chrM should map to GTF MT."""
        aliases = build_chrom_aliases(["chrM"], ["MT"])
        assert "MT" in aliases.get("chrM", [])

    def test_exact_match_kept(self):
        """Exact matching chrom names should be included."""
        aliases = build_chrom_aliases(["chr1"], ["chr1", "chr2"])
        assert "chr1" in aliases.get("chr1", [])


class TestBiotypeMap:

    def test_known_biotypes_mapped(self):
        assert BIOTYPE_MAP["miRNA"]      == "miRNA"
        assert BIOTYPE_MAP["pre_miRNA"]  == "miRNA"
        assert BIOTYPE_MAP["tRNA"]       == "tRNA"
        assert BIOTYPE_MAP["Mt_tRNA"]    == "tRNA"
        assert BIOTYPE_MAP["snoRNA"]     == "snoRNA"
        assert BIOTYPE_MAP["rRNA"]       == "rRNA"
        assert BIOTYPE_MAP["lncRNA"]     == "lncRNA"

    def test_priority_order(self):
        assert BIOTYPE_PRIORITY["miRNA"]  < BIOTYPE_PRIORITY["tRNA"]
        assert BIOTYPE_PRIORITY["tRNA"]   < BIOTYPE_PRIORITY["snoRNA"]
        assert BIOTYPE_PRIORITY["snoRNA"] < BIOTYPE_PRIORITY["rRNA"]
        assert BIOTYPE_PRIORITY["rRNA"]   < BIOTYPE_PRIORITY["lncRNA"]
        assert BIOTYPE_PRIORITY["lncRNA"] < BIOTYPE_PRIORITY["other"]
