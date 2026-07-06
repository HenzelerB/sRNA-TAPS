# -*- coding: utf-8 -*-
"""Tests for truth-set evaluation helpers."""

import pandas as pd
import pytest

from srnataps.evaluate import load_calls, run_evaluation


def test_run_evaluation_scores_truth_sites_with_coordinate_conversion(tmp_path):
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "chrom\tgenomic_pos\tstrand\tbiotype\tgene_id\tplanted_rate\n"
        "chr1\t101\t+\tmiRNA\tgeneA\t0.8\n"
        "chr1\t201\t+\tmiRNA\tgeneB\t0.7\n"
    )

    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\tcell_line\tfastq\n"
        "treat_HEK_R1\ttreat\tHEK\t/a.fq.gz\n"
        "pb_Ctrl_HEK_R1\tpb_Ctrl\tHEK\t/b.fq.gz\n"
    )

    calls_dir = tmp_path / "07.taps_calls" / "miRNA"
    calls_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            # True positive: truth genomic_pos 101 -> caller start 100.
            {
                "chrom": "chr1",
                "start": 100,
                "end": 101,
                "strand": "+",
                "context": "ACA",
                "mod_count": 8,
                "unmod_count": 2,
                "coverage": 10,
                "mod_rate": 0.8,
                "pvalue": 0.001,
                "padj": 0.01,
                "snp_flag": "PASS",
            },
            # False positive in treat.
            {
                "chrom": "chr1",
                "start": 300,
                "end": 301,
                "strand": "+",
                "context": "ACA",
                "mod_count": 8,
                "unmod_count": 2,
                "coverage": 10,
                "mod_rate": 0.8,
                "pvalue": 0.001,
                "padj": 0.01,
                "snp_flag": "PASS",
            },
        ]
    ).to_csv(calls_dir / "treat_HEK_R1_miRNA_taps.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            # This overlaps truth, but should not count when scoring treat only.
            {
                "chrom": "chr1",
                "start": 200,
                "end": 201,
                "strand": "+",
                "context": "ACA",
                "mod_count": 8,
                "unmod_count": 2,
                "coverage": 10,
                "mod_rate": 0.8,
                "pvalue": 0.001,
                "padj": 0.01,
                "snp_flag": "PASS",
            },
        ]
    ).to_csv(calls_dir / "pb_Ctrl_HEK_R1_miRNA_taps.tsv", sep="\t", index=False)

    summary, by_biotype = run_evaluation(
        truth_path=truth,
        calls_dir=tmp_path / "07.taps_calls",
        samples_tsv=samples,
        outdir=tmp_path / "10.truth_evaluation",
        max_padj=0.05,
        condition="treat",
    )

    row = summary.iloc[0]
    assert row["truth_sites"] == 2
    assert row["called_sites"] == 2
    assert row["true_positive_sites"] == 1
    assert row["false_positive_sites"] == 1
    assert row["false_negative_sites"] == 1
    assert row["precision"] == 0.5
    assert row["recall"] == 0.5

    bt = by_biotype.iloc[0]
    assert bt["biotype"] == "miRNA"
    assert bt["true_positive_sites"] == 1

    recovered = pd.read_csv(
        tmp_path / "10.truth_evaluation" / "truth_sites_recovery.tsv", sep="\t"
    )
    assert recovered["recovered"].tolist() == [True, False]


def test_load_calls_fails_when_calls_dir_is_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="calls directory not found"):
        load_calls(tmp_path / "missing_calls")


def test_overall_detection_is_independent_of_biotype_label(tmp_path):
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "chrom\tgenomic_pos\tstrand\tbiotype\tgene_id\tplanted_rate\n"
        "chr1\t101\t+\tsnoRNA\tgeneA\t0.8\n"
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\tcell_line\tfastq\n" "treat_HEK_R1\ttreat\tHEK\t/a.fq.gz\n"
    )
    calls_dir = tmp_path / "07.taps_calls" / "miRNA"
    calls_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 100,
                "end": 101,
                "strand": "+",
                "context": "ACA",
                "mod_count": 8,
                "unmod_count": 2,
                "coverage": 10,
                "mod_rate": 0.8,
                "pvalue": 0.001,
                "padj": 0.01,
                "snp_flag": "PASS",
            }
        ]
    ).to_csv(calls_dir / "treat_HEK_R1_miRNA_taps.tsv", sep="\t", index=False)

    summary, by_biotype = run_evaluation(
        truth_path=truth,
        calls_dir=tmp_path / "07.taps_calls",
        samples_tsv=samples,
        outdir=tmp_path / "evaluation",
        condition="treat",
    )

    assert summary.iloc[0]["recall"] == 1.0
    sno = by_biotype[by_biotype["biotype"] == "snoRNA"].iloc[0]
    assert sno["recall"] == 0.0
    recovery = pd.read_csv(
        tmp_path / "evaluation" / "truth_sites_recovery.tsv", sep="\t"
    )
    assert recovery.iloc[0]["recovered"]
    assert not recovery.iloc[0]["biotype_match"]
    assert recovery.iloc[0]["detected_biotypes"] == "miRNA"
