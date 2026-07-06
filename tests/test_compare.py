# -*- coding: utf-8 -*-
"""tests/test_compare.py - Tests for comparison grouping helpers."""

from srnataps.compare import discover_condition_prefixes


def test_discover_condition_prefixes_from_samples_tsv(tmp_path, monkeypatch):
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\tcell_line\tfastq\n"
        "untreated_A549_R1\tuntreated\tA549\t/a.fq.gz\n"
        "untreated_A549_R2\tuntreated\tA549\t/b.fq.gz\n"
        "pbOnly_A549_R1\tPB only\tA549\t/c.fq.gz\n"
        "pbOnly_A549_R2\tPB only\tA549\t/d.fq.gz\n"
        "tetpb_A549_R1\tTET+PB\tA549\t/e.fq.gz\n"
        "tetpb_A549_R2\tTET+PB\tA549\t/f.fq.gz\n"
    )
    monkeypatch.chdir(tmp_path)

    assert discover_condition_prefixes() == [
        "untreated_A549",
        "pbOnly_A549",
        "tetpb_A549",
    ]
