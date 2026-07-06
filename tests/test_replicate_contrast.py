from pathlib import Path

import pandas as pd

from srnataps.replicate_contrast import (
    beta_binomial_lrt,
    confidence_status,
    one_sided_mann_whitney,
    sample_from_bam_path,
    run_replicate_contrast,
)


def test_mann_whitney_detects_consistent_treatment_shift():
    treat = [0.7, 0.8, 0.75, 0.9, 0.85]
    control = [0.01, 0.02, 0.0, 0.03, 0.01]
    assert one_sided_mann_whitney(treat, control, 3) < 0.05


def test_mann_whitney_requires_minimum_replicates():
    assert one_sided_mann_whitney([0.8, 0.9], [0.0, 0.1], 3) == 1.0


def test_confidence_status_distinguishes_missing_from_negative_evidence():
    assert confidence_status(False, False, False) == "excluded_before_replicate_test"
    assert confidence_status(True, False, False) == "insufficient_replicate_evidence"
    assert confidence_status(True, True, False) == "replicate_not_supported"
    assert confidence_status(True, True, True) == "replicate_supported"


def test_sample_name_is_parsed_from_biotype_bam():
    path = Path("05.biotype_bams_3letter/miRNA/treat_HEK_R1_miRNA.sorted.bam")
    assert sample_from_bam_path(path, "miRNA") == "treat_HEK_R1"


def test_beta_binomial_detects_count_supported_shift():
    treat = [(70, 30), (80, 20), (75, 25), (85, 15), (78, 22)]
    control = [(2, 98), (1, 99), (3, 97), (2, 98), (1, 99)]
    assert beta_binomial_lrt(treat, control, 3) < 0.05


def test_beta_binomial_is_one_sided_and_requires_replicates():
    high = [(80, 20), (75, 25), (85, 15)]
    low = [(2, 98), (1, 99), (3, 97)]
    assert beta_binomial_lrt(low, high, 3) == 1.0
    assert beta_binomial_lrt(high[:2], low[:2], 3) == 1.0


def test_no_coverage_qualified_candidates_writes_empty_outputs(tmp_path):
    audit = pd.DataFrame(
        {
            "chrom": ["1"],
            "start": [10],
            "strand": ["+"],
            "coverage": [1],
            "pb_ctrl_coverage": [1],
            "no_treat_coverage": [1],
        }
    )
    audit_path = tmp_path / "audit.tsv"
    samples_path = tmp_path / "samples.tsv"
    out_path = tmp_path / "passing.tsv"
    all_out_path = tmp_path / "all.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    pd.DataFrame(
        {"sample": ["treat_R1"], "condition": ["treat"], "cell_line": ["HEK"]}
    ).to_csv(samples_path, sep="\t", index=False)

    _, passing = run_replicate_contrast(
        audit_path,
        tmp_path / "calls",
        samples_path,
        "miRNA",
        out_path,
        all_out_path,
    )

    assert passing.empty
    result = pd.read_csv(all_out_path, sep="\t")
    assert result.loc[0, "confidence_status"] == "excluded_before_replicate_test"
