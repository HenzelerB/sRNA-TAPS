import pandas as pd

from srnataps.tier import filter_stringent_calls


def test_stringent_filter_requires_depth_effect_and_both_qvalues(tmp_path):
    audit = pd.DataFrame(
        {
            "coverage": [100, 100, 2],
            "pb_ctrl_coverage": [100, 100, 100],
            "no_treat_coverage": [100, 100, 100],
            "replicate_delta_pb_ctrl": [0.5, 0.5, 0.5],
            "replicate_delta_no_treat": [0.5, 0.5, 0.5],
            "replicate_pb_padj": [1e-25, 1e-25, 1e-25],
            "replicate_no_treat_padj": [1e-25, 0.1, 1e-25],
            "pb_available": [True, True, True],
            "no_treat_available": [True, True, True],
        }
    )
    audit_path = tmp_path / "audit.tsv"
    out_path = tmp_path / "stringent.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)

    passing = filter_stringent_calls(audit_path, out_path)

    assert len(passing) == 1
    assert passing.iloc[0]["evidence_tier"] == "stringent"
