import pandas as pd

from srnataps.contrast import contrast_pooled_calls


def calls(rows):
    return pd.DataFrame(rows)


def test_control_contrast_keeps_treatment_specific_signal():
    treat = calls(
        [
            {
                "chrom": "1",
                "start": 10,
                "end": 11,
                "strand": "+",
                "mod_count": 80,
                "unmod_count": 20,
                "coverage": 100,
                "mod_rate": 0.8,
            },
            {
                "chrom": "1",
                "start": 20,
                "end": 21,
                "strand": "+",
                "mod_count": 20,
                "unmod_count": 80,
                "coverage": 100,
                "mod_rate": 0.2,
            },
        ]
    )
    pb = calls(
        [
            {
                "chrom": "1",
                "start": 10,
                "strand": "+",
                "pb_ctrl_mod_count": 2,
                "pb_ctrl_unmod_count": 98,
                "pb_ctrl_coverage": 100,
                "pb_ctrl_mod_rate": 0.02,
            },
            {
                "chrom": "1",
                "start": 20,
                "strand": "+",
                "pb_ctrl_mod_count": 18,
                "pb_ctrl_unmod_count": 82,
                "pb_ctrl_coverage": 100,
                "pb_ctrl_mod_rate": 0.18,
            },
        ]
    )
    untreated = calls(
        [
            {
                "chrom": "1",
                "start": 10,
                "strand": "+",
                "no_treat_mod_count": 1,
                "no_treat_unmod_count": 99,
                "no_treat_coverage": 100,
                "no_treat_mod_rate": 0.01,
            },
            {
                "chrom": "1",
                "start": 20,
                "strand": "+",
                "no_treat_mod_count": 19,
                "no_treat_unmod_count": 81,
                "no_treat_coverage": 100,
                "no_treat_mod_rate": 0.19,
            },
        ]
    )

    result = contrast_pooled_calls(treat, pb, untreated)

    assert bool(result.loc[result["start"] == 10, "contrast_pass"].iloc[0])
    assert not bool(result.loc[result["start"] == 20, "contrast_pass"].iloc[0])


def test_control_contrast_requires_control_coverage():
    treat = calls(
        [
            {
                "chrom": "1",
                "start": 10,
                "end": 11,
                "strand": "+",
                "mod_count": 10,
                "unmod_count": 0,
                "coverage": 10,
                "mod_rate": 1.0,
            }
        ]
    )
    empty_pb = calls(
        [
            {
                "chrom": "2",
                "start": 20,
                "strand": "+",
                "pb_ctrl_mod_count": 0,
                "pb_ctrl_unmod_count": 10,
                "pb_ctrl_coverage": 10,
                "pb_ctrl_mod_rate": 0.0,
            }
        ]
    )
    empty_untreated = calls(
        [
            {
                "chrom": "2",
                "start": 20,
                "strand": "+",
                "no_treat_mod_count": 0,
                "no_treat_unmod_count": 10,
                "no_treat_coverage": 10,
                "no_treat_mod_rate": 0.0,
            }
        ]
    )

    result = contrast_pooled_calls(treat, empty_pb, empty_untreated)

    assert not bool(result.iloc[0]["contrast_pass"])
