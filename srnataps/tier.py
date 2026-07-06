"""Filter replicate evidence audits into a configurable stringent tier."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def filter_stringent_calls(
    audit_path,
    out_path,
    maximum_padj=1e-20,
    minimum_delta=0.1,
    minimum_coverage=5,
):
    calls = pd.read_csv(audit_path, sep="\t", low_memory=False)
    if calls.empty:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        calls.to_csv(out_path, sep="\t", index=False)
        return calls
    required = {
        "coverage",
        "pb_ctrl_coverage",
        "no_treat_coverage",
        "replicate_delta_pb_ctrl",
        "replicate_delta_no_treat",
        "replicate_pb_padj",
        "replicate_no_treat_padj",
        "pb_available",
        "no_treat_available",
    }
    missing = required - set(calls.columns)
    if missing:
        raise ValueError(f"{audit_path} missing columns: {', '.join(sorted(missing))}")

    passing = calls[
        calls["pb_available"].fillna(False)
        & calls["no_treat_available"].fillna(False)
        & (calls["coverage"] >= minimum_coverage)
        & (calls["pb_ctrl_coverage"] >= minimum_coverage)
        & (calls["no_treat_coverage"] >= minimum_coverage)
        & (calls["replicate_delta_pb_ctrl"] >= minimum_delta)
        & (calls["replicate_delta_no_treat"] >= minimum_delta)
        & (calls["replicate_pb_padj"] <= maximum_padj)
        & (calls["replicate_no_treat_padj"] <= maximum_padj)
    ].copy()
    passing["evidence_tier"] = "stringent"

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    passing.to_csv(out_path, sep="\t", index=False)
    return passing


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--maximum-padj", type=float, default=1e-20)
    parser.add_argument("--minimum-delta", type=float, default=0.1)
    parser.add_argument("--minimum-coverage", type=float, default=5)
    return parser.parse_args()


def main():
    args = parse_args()
    passing = filter_stringent_calls(
        args.audit,
        args.out,
        args.maximum_padj,
        args.minimum_delta,
        args.minimum_coverage,
    )
    print(f"Stringent sites: {len(passing):,}")


if __name__ == "__main__":
    main()
