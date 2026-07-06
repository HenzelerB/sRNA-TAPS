# -*- coding: utf-8 -*-
"""Control-aware testing for pooled sRNA-TAPS condition calls."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

REQUIRED = {
    "chrom",
    "start",
    "end",
    "strand",
    "mod_count",
    "unmod_count",
    "coverage",
    "mod_rate",
}
KEY = ["chrom", "start", "strand"]


def load_calls(path, prefix=None):
    calls = pd.read_csv(path, sep="\t")
    missing = REQUIRED - set(calls.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")

    calls = calls.copy()
    calls["chrom"] = calls["chrom"].astype(str)
    calls["start"] = pd.to_numeric(calls["start"], errors="raise").astype(int)
    for column in ["mod_count", "unmod_count", "coverage", "mod_rate"]:
        calls[column] = pd.to_numeric(calls[column], errors="raise")

    if calls.duplicated(KEY).any():
        raise ValueError(f"{path} contains duplicate genomic site keys")
    if prefix is None:
        return calls

    return calls[KEY + ["mod_count", "unmod_count", "coverage", "mod_rate"]].rename(
        columns={
            "mod_count": f"{prefix}_mod_count",
            "unmod_count": f"{prefix}_unmod_count",
            "coverage": f"{prefix}_coverage",
            "mod_rate": f"{prefix}_mod_rate",
        }
    )


def _one_sided_fisher(treat_mod, treat_unmod, control_mod, control_unmod):
    table = [
        [int(round(treat_mod)), int(round(treat_unmod))],
        [int(round(control_mod)), int(round(control_unmod))],
    ]
    if sum(table[0]) == 0 or sum(table[1]) == 0:
        return 1.0
    return float(fisher_exact(table, alternative="greater").pvalue)


def _adjust(pvalues):
    if not pvalues:
        return []
    return multipletests(pvalues, method="fdr_bh")[1]


def contrast_pooled_calls(
    treat,
    pb_control,
    no_treat,
    min_treat_coverage=5,
    min_control_coverage=5,
    min_delta=0.1,
    max_padj=0.05,
):
    """Compare pooled treatment counts with matched pooled controls."""
    result = treat.merge(pb_control, on=KEY, how="left", validate="one_to_one")
    result = result.merge(no_treat, on=KEY, how="left", validate="one_to_one")

    control_columns = [
        f"{prefix}_{field}"
        for prefix in ("pb_ctrl", "no_treat")
        for field in ("mod_count", "unmod_count", "coverage", "mod_rate")
    ]
    result[control_columns] = result[control_columns].fillna(0.0)

    result["delta_pb_ctrl"] = result["mod_rate"] - result["pb_ctrl_mod_rate"]
    result["delta_no_treat"] = result["mod_rate"] - result["no_treat_mod_rate"]

    pb_pvalues = []
    untreated_pvalues = []
    for row in result.itertuples(index=False):
        pb_pvalues.append(
            _one_sided_fisher(
                row.mod_count,
                row.unmod_count,
                row.pb_ctrl_mod_count,
                row.pb_ctrl_unmod_count,
            )
        )
        untreated_pvalues.append(
            _one_sided_fisher(
                row.mod_count,
                row.unmod_count,
                row.no_treat_mod_count,
                row.no_treat_unmod_count,
            )
        )

    result["pvalue"] = pb_pvalues
    result["padj"] = _adjust(pb_pvalues)
    result["no_treat_pvalue"] = untreated_pvalues
    result["no_treat_padj"] = _adjust(untreated_pvalues)

    result["contrast_pass"] = (
        (result["coverage"] >= min_treat_coverage)
        & (result["pb_ctrl_coverage"] >= min_control_coverage)
        & (result["no_treat_coverage"] >= min_control_coverage)
        & (result["delta_pb_ctrl"] >= min_delta)
        & (result["delta_no_treat"] >= min_delta)
        & (result["padj"] <= max_padj)
        & (result["no_treat_padj"] <= max_padj)
    )
    return result


def run_contrast(
    treat_path,
    pb_control_path,
    no_treat_path,
    out_path,
    all_out_path=None,
    **kwargs,
):
    treat = load_calls(treat_path)
    pb_control = load_calls(pb_control_path, "pb_ctrl")
    no_treat = load_calls(no_treat_path, "no_treat")
    result = contrast_pooled_calls(treat, pb_control, no_treat, **kwargs)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    passing = result[result["contrast_pass"]].copy()
    passing.to_csv(out_path, sep="\t", index=False)

    if all_out_path:
        all_out_path = Path(all_out_path)
        all_out_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(all_out_path, sep="\t", index=False)
    return result, passing


def parse_args():
    parser = argparse.ArgumentParser(
        description="Contrast pooled treat calls with pb_ctrl and no_treat"
    )
    parser.add_argument("--treat", required=True)
    parser.add_argument("--pb-control", required=True)
    parser.add_argument("--no-treat", required=True)
    parser.add_argument("--out", required=True, help="passing contrast calls")
    parser.add_argument("--all-out", help="optional table containing all tested sites")
    parser.add_argument("--min-treat-coverage", type=float, default=5)
    parser.add_argument("--min-control-coverage", type=float, default=5)
    parser.add_argument("--min-delta", type=float, default=0.1)
    parser.add_argument("--max-padj", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    result, passing = run_contrast(
        treat_path=args.treat,
        pb_control_path=args.pb_control,
        no_treat_path=args.no_treat,
        out_path=args.out,
        all_out_path=args.all_out,
        min_treat_coverage=args.min_treat_coverage,
        min_control_coverage=args.min_control_coverage,
        min_delta=args.min_delta,
        max_padj=args.max_padj,
    )
    print(f"Tested sites: {len(result):,}")
    print(f"Passing sites: {len(passing):,}")


if __name__ == "__main__":
    main()
