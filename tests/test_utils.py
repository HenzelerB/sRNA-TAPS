# -*- coding: utf-8 -*-
"""tests/test_utils.py — Tests for srnataps.utils"""

import pytest
from srnataps.utils import (
    detect_cell_line,
    detect_condition,
    get_chroms,
    is_condition,
    normalize_condition,
)


class TestDetectCellLine:
    def test_hek(self):       assert detect_cell_line("treat_HEK_R1")        == "HEK"
    def test_caco2(self):     assert detect_cell_line("no-treat_Ctrl_Caco2") == "Caco2"
    def test_unknown(self):   assert detect_cell_line("some_other_sample")   == "unknown"


class TestDetectCondition:
    def test_notreat(self):   assert detect_condition("no-treat_Ctrl_HEK_R1")  == "no_treat"
    def test_pb(self):        assert detect_condition("pb_Ctrl_Caco2_R2")       == "pb_ctrl"
    def test_treat(self):     assert detect_condition("treat_HEK_R3")           == "treat"
    def test_unknown(self):   assert detect_condition("mystery_sample")         == "unknown"


class TestNormalizeCondition:
    @pytest.mark.parametrize("label", [
        "no_treat", "no-treat", "notreat", "No Treat", "untreated", "untr",
    ])
    def test_notreat_aliases(self, label):
        assert normalize_condition(label) == "no_treat"
        assert is_condition(label, "no_treat")

    @pytest.mark.parametrize("label", [
        "pb_ctrl", "pb_Ctrl", "PB-Ctrl", "pb control", "PB only", "pb",
    ])
    def test_pb_aliases(self, label):
        assert normalize_condition(label) == "pb_ctrl"
        assert is_condition(label, "pb_ctrl")

    @pytest.mark.parametrize("label", [
        "treat", "treated", "TET+PB", "tet_pb", "full_taps",
    ])
    def test_treat_aliases(self, label):
        assert normalize_condition(label) == "treat"
        assert is_condition(label, "treat")

    def test_unknown_passthrough(self):
        assert normalize_condition("custom_condition") == "custom_condition"

    def test_unknown_strict_raises(self):
        with pytest.raises(ValueError):
            normalize_condition("custom_condition", strict=True)


class TestGetChroms:
    def test_returns_chroms(self, bam_notreat):
        chroms = get_chroms(bam_notreat)
        assert isinstance(chroms, list)
        assert len(chroms) >= 1
        assert "chr1" in chroms
