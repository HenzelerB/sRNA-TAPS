"""Summarize condition-level sRNA-TAPS evidence tiers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def count_rows(path):
    try:
        return len(pd.read_csv(path, sep="\t", low_memory=False))
    except pd.errors.EmptyDataError:
        return 0


def build_summary(cell_lines, biotypes, directories, out_path):
    templates = {
        "pooled_test_universe": "{biotype}/treat_{cell_line}_pooled_{biotype}_taps.tsv",
        "control_contrast": "{biotype}/treat_{cell_line}_contrast_{biotype}_taps.tsv",
        "replicate_discovery": "{biotype}/treat_{cell_line}_replicate_{biotype}_taps.tsv",
        "stringent": "{biotype}/treat_{cell_line}_stringent_{biotype}_taps.tsv",
    }
    rows = []
    for cell_line in cell_lines:
        for biotype in biotypes:
            for tier, directory in directories.items():
                path = Path(directory) / templates[tier].format(
                    cell_line=cell_line, biotype=biotype
                )
                rows.append(
                    {
                        "cell_line": cell_line,
                        "biotype": biotype,
                        "evidence_tier": tier,
                        "sites": count_rows(path),
                        "path": str(path),
                    }
                )
    summary = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, sep="\t", index=False)
    return summary


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-lines", nargs="+", required=True)
    parser.add_argument("--biotypes", nargs="+", required=True)
    parser.add_argument("--pooled-dir", required=True)
    parser.add_argument("--contrast-dir", required=True)
    parser.add_argument("--replicate-dir", required=True)
    parser.add_argument("--stringent-dir", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    summary = build_summary(
        args.cell_lines,
        args.biotypes,
        {
            "pooled_test_universe": args.pooled_dir,
            "control_contrast": args.contrast_dir,
            "replicate_discovery": args.replicate_dir,
            "stringent": args.stringent_dir,
        },
        args.out,
    )
    print(f"Evidence summary rows: {len(summary):,}")


if __name__ == "__main__":
    main()
