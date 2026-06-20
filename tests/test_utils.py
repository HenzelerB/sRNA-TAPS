# -*- coding: utf-8 -*-
"""tests/test_utils.py — Tests for srnataps.utils"""

import pytest
from srnataps.utils import detect_cell_line, detect_condition, get_chroms


class TestDetectCellLine:
    def test_hek(self):       assert detect_cell_line("treat_HEK_R1")        == "HEK"
    def test_caco2(self):     assert detect_cell_line("no-treat_Ctrl_Caco2") == "Caco2"
    def test_unknown(self):   assert detect_cell_line("some_other_sample")   == "unknown"


class TestDetectCondition:
    def test_notreat(self):   assert detect_condition("no-treat_Ctrl_HEK_R1")  == "no_treat"
    def test_pb(self):        assert detect_condition("pb_Ctrl_Caco2_R2")       == "pb_ctrl"
    def test_treat(self):     assert detect_condition("treat_HEK_R3")           == "treat"
    def test_unknown(self):   assert detect_condition("mystery_sample")         == "unknown"


class TestGetChroms:
    def test_returns_chroms(self, bam_notreat):
        chroms = get_chroms(bam_notreat)
        assert isinstance(chroms, list)
        assert len(chroms) >= 1
        assert "chr1" in chroms
