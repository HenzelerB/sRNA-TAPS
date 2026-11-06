from pathlib import Path

import pandas as pd

from srnataps.replicate_contrast import (
    bayesian_shift_probability,
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


def test_bayesian_shift_probability_detects_supported_delta():
    treat = [(70, 30), (80, 20), (75, 25), (85, 15)]
    control = [(2, 98), (1, 99), (3, 97), (2, 98)]
    assert bayesian_shift_probability(treat, control, 3, minimum_delta=0.1) > 0.99


def test_bayesian_shift_probability_is_directional_and_requires_replicates():
    high = [(80, 20), (75, 25), (85, 15)]
    low = [(2, 98), (1, 99), (3, 97)]
    assert bayesian_shift_probability(low, high, 3, minimum_delta=0.1) < 0.01
    assert bayesian_shift_probability(high[:2], low[:2], 3, minimum_delta=0.1) == 0.0


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


def test_bayesian_no_treat_primary_keeps_pb_as_reported_background(tmp_path):
    audit = pd.DataFrame(
        {
            "chrom": ["1"],
            "start": [10],
            "strand": ["+"],
            "coverage": [300],
            "pb_ctrl_coverage": [300],
            "no_treat_coverage": [300],
            "mod_rate": [0.8],
            "pb_ctrl_mod_rate": [0.75],
            "no_treat_mod_rate": [0.01],
            "delta_pb_ctrl": [0.05],
            "delta_no_treat": [0.79],
            "pvalue": [0.001],
            "padj": [0.001],
        }
    )
    audit_path = tmp_path / "audit.tsv"
    samples_path = tmp_path / "samples.tsv"
    calls_dir = tmp_path / "calls" / "miRNA"
    out_path = tmp_path / "passing.tsv"
    all_out_path = tmp_path / "all.tsv"
    calls_dir.mkdir(parents=True)
    audit.to_csv(audit_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "sample": [
                "treat_R1",
                "treat_R2",
                "treat_R3",
                "pb_R1",
                "pb_R2",
                "pb_R3",
                "no_R1",
                "no_R2",
                "no_R3",
            ],
            "condition": [
                "treat",
                "treat",
                "treat",
                "pb_ctrl",
                "pb_ctrl",
                "pb_ctrl",
                "no_treat",
                "no_treat",
                "no_treat",
            ],
            "cell_line": ["HEK"] * 9,
        }
    ).to_csv(samples_path, sep="\t", index=False)

    for sample, modified, unmodified in [
        ("treat_R1", 80, 20),
        ("treat_R2", 81, 19),
        ("treat_R3", 79, 21),
        ("pb_R1", 75, 25),
        ("pb_R2", 74, 26),
        ("pb_R3", 76, 24),
        ("no_R1", 1, 99),
        ("no_R2", 2, 98),
        ("no_R3", 1, 99),
    ]:
        pd.DataFrame(
            {
                "chrom": ["1"],
                "start": [10],
                "end": [11],
                "strand": ["+"],
                "mod_count": [modified],
                "unmod_count": [unmodified],
                "coverage": [modified + unmodified],
                "mod_rate": [modified / (modified + unmodified)],
            }
        ).to_csv(calls_dir / f"{sample}_miRNA_taps.tsv", sep="\t", index=False)

    _, passing = run_replicate_contrast(
        audit_path,
        tmp_path / "calls",
        samples_path,
        "miRNA",
        out_path,
        all_out_path,
        statistical_test="bayesian",
        primary_control="no_treat",
        minimum_posterior_probability=0.95,
    )

    assert len(passing) == 1
    assert passing.iloc[0]["primary_control"] == "no_treat"
    assert passing.iloc[0]["replicate_pb_ctrl_posterior_probability"] < 0.95
    assert passing.iloc[0]["replicate_no_treat_posterior_probability"] > 0.95
