# -*- coding: utf-8 -*-
"""
Truth-set evaluation for simulated sRNA-TAPS data.

The simulator writes truth.tsv with 1-based genomic positions. The caller writes
BED-like 0-based start / 1-based end coordinates. This module keeps that
conversion explicit and scores high-confidence calls against planted sites.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from srnataps.utils import detect_condition, normalize_condition

TRUTH_REQUIRED = {
    "chrom",
    "genomic_pos",
    "strand",
    "biotype",
    "gene_id",
    "planted_rate",
}
CALL_REQUIRED = {"chrom", "start", "strand", "mod_rate", "coverage", "padj"}


def site_key(chrom, start0, strand, biotype):
    return (str(chrom), int(start0), str(strand), str(biotype))


def detection_key(chrom, start0, strand):
    """Genomic modification-site identity, independent of annotation label."""
    return (str(chrom), int(start0), str(strand))


def load_sample_conditions(samples_tsv):
    """Return sample -> canonical condition from samples.tsv, or an empty dict."""
    if not samples_tsv:
        return {}
    path = Path(samples_tsv)
    if not path.exists():
        return {}

    samples = pd.read_csv(path, sep="\t")
    if not {"sample", "condition"}.issubset(samples.columns):
        return {}

    out = {}
    for _, row in samples.iterrows():
        out[str(row["sample"])] = normalize_condition(row["condition"], strict=False)
    return out


def load_truth(path):
    truth = pd.read_csv(path, sep="\t", dtype={"chrom": "string"}, low_memory=False)
    missing = TRUTH_REQUIRED - set(truth.columns)
    if missing:
        raise ValueError(f"truth file missing columns: {', '.join(sorted(missing))}")

    truth = truth.copy()
    truth["start"] = truth["genomic_pos"].astype(int) - 1
    truth["truth_key"] = [
        site_key(r.chrom, r.start, r.strand, r.biotype) for r in truth.itertuples()
    ]
    truth["detection_key"] = [
        detection_key(r.chrom, r.start, r.strand) for r in truth.itertuples()
    ]
    return truth


def parse_call_filename(path):
    """
    Parse <sample>_<biotype>_taps.tsv or <sample>_<biotype>_taps_annotated.tsv.

    Sample names and biotypes can contain underscores, so this uses the parent
    directory as the biotype and removes that exact suffix from the filename.
    """
    biotype = path.parent.name
    patterns = [
        rf"^(?P<sample>.+)_{re.escape(biotype)}_taps(?:_annotated)?\.tsv$",
        r"^(?P<sample>.+)_taps(?:_annotated)?\.tsv$",
    ]
    for pattern in patterns:
        match = re.match(pattern, path.name)
        if match:
            return match.group("sample"), biotype
    return path.stem, biotype


def iter_call_files(calls_dir, annotated=False):
    calls_dir = Path(calls_dir)
    if not calls_dir.exists():
        raise FileNotFoundError(
            f"calls directory not found: {calls_dir}. "
            "Run the TAPS calling stage first, or pass the correct --calls-dir."
        )
    suffix = "_taps_annotated.tsv" if annotated else "_taps.tsv"
    return sorted(calls_dir.glob(f"*/*{suffix}"))


def load_calls(calls_dir, samples_tsv=None, annotated=False):
    sample_conditions = load_sample_conditions(samples_tsv)
    frames = []
    call_files = iter_call_files(calls_dir, annotated=annotated)
    if not call_files:
        suffix = "_taps_annotated.tsv" if annotated else "_taps.tsv"
        raise ValueError(f"no {suffix} files found below calls directory: {calls_dir}")

    for path in call_files:
        sample, biotype = parse_call_filename(path)
        try:
            df = pd.read_csv(
                path, sep="\t", dtype={"chrom": "string"}, low_memory=False
            )
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue

        missing = CALL_REQUIRED - set(df.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")

        df = df.copy()
        df["sample"] = sample
        df["biotype"] = biotype
        df["condition"] = (
            sample_conditions.get(sample) or detect_condition(sample) or "unknown"
        )
        df["condition"] = df["condition"].map(
            lambda v: normalize_condition(v, strict=False)
        )
        df["source_file"] = str(path)
        df["call_key"] = [
            site_key(r.chrom, r.start, r.strand, r.biotype) for r in df.itertuples()
        ]
        df["detection_key"] = [
            detection_key(r.chrom, r.start, r.strand) for r in df.itertuples()
        ]
        frames.append(df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def filter_calls(
    calls, min_coverage=1, min_mod_rate=0.0, max_padj=1.0, condition="treat"
):
    if calls.empty:
        return calls

    out = calls.copy()
    out["coverage"] = pd.to_numeric(out["coverage"], errors="coerce")
    out["mod_rate"] = pd.to_numeric(out["mod_rate"], errors="coerce")
    out["padj"] = pd.to_numeric(out["padj"], errors="coerce")

    if condition != "all":
        out = out[out["condition"] == normalize_condition(condition, strict=False)]

    return out[
        (out["coverage"] >= min_coverage)
        & (out["mod_rate"] >= min_mod_rate)
        & (out["padj"] <= max_padj)
    ].copy()


def summarise_site_level(
    truth,
    calls,
    truth_key_column="detection_key",
    call_key_column="detection_key",
):
    truth_keys = set(truth[truth_key_column])
    call_keys = set(calls[call_key_column]) if not calls.empty else set()

    tp = len(truth_keys & call_keys)
    fp = len(call_keys - truth_keys)
    fn = len(truth_keys - call_keys)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0

    return {
        "truth_sites": len(truth_keys),
        "called_sites": len(call_keys),
        "true_positive_sites": tp,
        "false_positive_sites": fp,
        "false_negative_sites": fn,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
    }


def summarise_by_biotype(truth, calls):
    rows = []
    biotypes = sorted(
        set(truth["biotype"]) | set(calls["biotype"] if not calls.empty else [])
    )
    for biotype in biotypes:
        t = truth[truth["biotype"] == biotype]
        c = calls[calls["biotype"] == biotype] if not calls.empty else calls
        rows.append(
            {
                "biotype": biotype,
                **summarise_site_level(
                    t,
                    c,
                    truth_key_column="truth_key",
                    call_key_column="call_key",
                ),
            }
        )
    return pd.DataFrame(rows)


def build_truth_site_table(truth, calls):
    detection_keys = set(calls["detection_key"]) if not calls.empty else set()
    classification_keys = set(calls["call_key"]) if not calls.empty else set()
    detected_biotypes = {}
    if not calls.empty:
        for key, group in calls.groupby("detection_key"):
            detected_biotypes[key] = ",".join(sorted(set(group["biotype"])))

    out = truth.copy()
    out["recovered"] = out["detection_key"].map(lambda key: key in detection_keys)
    out["biotype_match"] = out["truth_key"].map(lambda key: key in classification_keys)
    out["detected_biotypes"] = out["detection_key"].map(
        lambda key: detected_biotypes.get(key, "")
    )
    return out.drop(columns=["truth_key", "detection_key"])


def build_call_site_table(calls, truth):
    if calls.empty:
        return calls
    truth_detection_keys = set(truth["detection_key"])
    truth_classification_keys = set(truth["truth_key"])
    out = calls.copy()
    out["truth_match"] = out["detection_key"].map(
        lambda key: key in truth_detection_keys
    )
    out["truth_biotype_match"] = out["call_key"].map(
        lambda key: key in truth_classification_keys
    )
    return out.drop(columns=["call_key", "detection_key"])


def run_evaluation(
    truth_path,
    calls_dir,
    outdir,
    samples_tsv=None,
    annotated=False,
    min_coverage=1,
    min_mod_rate=0.0,
    max_padj=1.0,
    condition="treat",
):
    truth = load_truth(truth_path)
    calls = load_calls(calls_dir, samples_tsv=samples_tsv, annotated=annotated)
    selected_calls = filter_calls(
        calls,
        min_coverage=min_coverage,
        min_mod_rate=min_mod_rate,
        max_padj=max_padj,
        condition=condition,
    )

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame([summarise_site_level(truth, selected_calls)])
    by_biotype = summarise_by_biotype(truth, selected_calls)
    truth_sites = build_truth_site_table(truth, selected_calls)
    call_sites = build_call_site_table(selected_calls, truth)

    summary.to_csv(outdir / "truth_evaluation_summary.tsv", sep="\t", index=False)
    by_biotype.to_csv(outdir / "truth_evaluation_by_biotype.tsv", sep="\t", index=False)
    truth_sites.to_csv(outdir / "truth_sites_recovery.tsv", sep="\t", index=False)
    call_sites.to_csv(outdir / "called_sites_truth_overlap.tsv", sep="\t", index=False)

    return summary, by_biotype


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate sRNA-TAPS calls against truth.tsv"
    )
    parser.add_argument("--truth", required=True, help="truth.tsv from the simulator")
    parser.add_argument(
        "--calls-dir", default="07.taps_calls", help="directory containing call TSVs"
    )
    parser.add_argument(
        "--samples-tsv", default="samples.tsv", help="sample sheet for condition labels"
    )
    parser.add_argument(
        "--outdir", default="10.truth_evaluation", help="output directory"
    )
    parser.add_argument(
        "--annotated", action="store_true", help="read *_taps_annotated.tsv files"
    )
    parser.add_argument("--min-coverage", type=float, default=1)
    parser.add_argument("--min-mod-rate", type=float, default=0.0)
    parser.add_argument("--max-padj", type=float, default=1.0)
    parser.add_argument(
        "--condition", default="treat", help="condition to score, or 'all'"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    summary, by_biotype = run_evaluation(
        truth_path=args.truth,
        calls_dir=args.calls_dir,
        samples_tsv=args.samples_tsv,
        outdir=args.outdir,
        annotated=args.annotated,
        min_coverage=args.min_coverage,
        min_mod_rate=args.min_mod_rate,
        max_padj=args.max_padj,
        condition=args.condition,
    )
    print(summary.to_string(index=False))
    print()
    print(by_biotype.to_string(index=False))


if __name__ == "__main__":
    main()
