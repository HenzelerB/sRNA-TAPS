import pandas as pd

from srnataps.evidence_summary import build_summary


def test_evidence_summary_counts_empty_and_populated_tiers(tmp_path):
    directories = {}
    filenames = {
        "pooled_test_universe": "treat_HEK_pooled_miRNA_taps.tsv",
        "control_contrast": "treat_HEK_contrast_miRNA_taps.tsv",
        "replicate_discovery": "treat_HEK_replicate_miRNA_taps.tsv",
        "stringent": "treat_HEK_stringent_miRNA_taps.tsv",
    }
    for tier, filename in filenames.items():
        directory = tmp_path / tier
        (directory / "miRNA").mkdir(parents=True)
        pd.DataFrame(
            {"start": [1, 2] if tier == "pooled_test_universe" else []}
        ).to_csv(directory / "miRNA" / filename, sep="\t", index=False)
        directories[tier] = directory

    summary = build_summary(["HEK"], ["miRNA"], directories, tmp_path / "summary.tsv")

    counts = dict(zip(summary.evidence_tier, summary.sites))
    assert counts["pooled_test_universe"] == 2
    assert counts["stringent"] == 0
