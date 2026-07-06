# -*- coding: utf-8 -*-
"""Replicate-aware condition testing for pooled sRNA-TAPS candidates."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, gammaln
from scipy.stats import chi2, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from srnataps.utils import normalize_condition

KEY = ["chrom", "start", "strand"]
CONFIDENCE_STATUSES = {
    "supported": "replicate_supported",
    "insufficient": "insufficient_replicate_evidence",
    "not_supported": "replicate_not_supported",
    "excluded": "excluded_before_replicate_test",
}


def load_samples(path, cell_line=None):
    samples = pd.read_csv(path, sep="\t", dtype={"sample": str})
    required = {"sample", "condition"}
    if not required.issubset(samples.columns):
        raise ValueError(f"{path} must contain sample and condition columns")
    if cell_line is not None:
        if "cell_line" not in samples.columns:
            raise ValueError(f"{path} must contain cell_line when --cell-line is used")
        samples = samples[samples["cell_line"].astype(str) == str(cell_line)]
    return {
        str(row.sample): normalize_condition(row.condition, strict=False)
        for row in samples.itertuples()
    }


def sample_from_call_path(path, biotype):
    suffix = f"_{biotype}_taps.tsv"
    if not path.name.endswith(suffix):
        raise ValueError(f"unexpected call filename: {path}")
    return path.name[: -len(suffix)]


def load_replicate_rates(calls_dir, biotype, sample_conditions, candidate_keys):
    rates = defaultdict(lambda: defaultdict(list))
    calls_dir = Path(calls_dir) / biotype

    for path in sorted(calls_dir.glob(f"*_{biotype}_taps.tsv")):
        sample = sample_from_call_path(path, biotype)
        condition = sample_conditions.get(sample)
        if condition not in {"treat", "pb_ctrl", "no_treat"}:
            continue
        try:
            calls = pd.read_csv(
                path, sep="\t", dtype={"chrom": "string"}, low_memory=False
            )
        except pd.errors.EmptyDataError:
            continue
        if calls.empty:
            continue

        calls["start"] = pd.to_numeric(calls["start"], errors="raise").astype(int)
        calls["mod_rate"] = pd.to_numeric(calls["mod_rate"], errors="raise")
        for row in calls.itertuples():
            key = (str(row.chrom), int(row.start), str(row.strand))
            if key in candidate_keys:
                if hasattr(row, "mod_count") and hasattr(row, "unmod_count"):
                    rates[key][condition].append(
                        (float(row.mod_count), float(row.unmod_count))
                    )
                else:
                    rate = float(row.mod_rate)
                    rates[key][condition].append((rate, 1.0 - rate))
    return rates


def sample_from_bam_path(path, biotype):
    suffix = f"_{biotype}.sorted.bam"
    if not path.name.endswith(suffix):
        raise ValueError(f"unexpected BAM filename: {path}")
    return path.name[: -len(suffix)]


def load_replicate_rates_from_bams(
    bam_dir,
    biotype,
    sample_conditions,
    candidate_keys,
    minimum_base_quality=20,
    minimum_mapping_quality=10,
    minimum_sample_coverage=1,
):
    """Collect candidate-site rates directly from each restored-sequence BAM."""
    import pysam

    rates = defaultdict(lambda: defaultdict(list))
    candidate_keys = set(candidate_keys)
    bam_dir = Path(bam_dir) / biotype

    for path in sorted(bam_dir.glob(f"*_{biotype}.sorted.bam")):
        sample = sample_from_bam_path(path, biotype)
        condition = sample_conditions.get(sample)
        if condition not in {"treat", "pb_ctrl", "no_treat"}:
            continue

        counts = defaultdict(lambda: [0.0, 0.0])
        with pysam.AlignmentFile(path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if (
                    read.is_unmapped
                    or read.query_sequence is None
                    or read.is_secondary
                    or read.is_supplementary
                    or read.mapping_quality < minimum_mapping_quality
                ):
                    continue

                try:
                    nh = read.get_tag("NH")
                except KeyError:
                    try:
                        nh = read.get_tag("XA")
                    except KeyError:
                        nh = 1
                weight = 1.0 / max(nh, 1)
                strand = "-" if read.is_reverse else "+"
                chrom = bam.get_reference_name(read.reference_id)
                qualities = read.query_qualities

                for query_pos, reference_pos in read.get_aligned_pairs(
                    matches_only=True
                ):
                    key = (str(chrom), int(reference_pos), strand)
                    if key not in candidate_keys:
                        continue
                    if (
                        qualities is not None
                        and qualities[query_pos] < minimum_base_quality
                    ):
                        continue

                    base = read.query_sequence[query_pos].upper()
                    unmodified, modified = ("C", "T") if strand == "+" else ("G", "A")
                    if base == modified:
                        counts[key][0] += weight
                    elif base == unmodified:
                        counts[key][1] += weight

        for key, (modified_count, unmodified_count) in counts.items():
            coverage = modified_count + unmodified_count
            if coverage >= minimum_sample_coverage:
                rates[key][condition].append((modified_count, unmodified_count))
    return rates


def observation_rates(observations):
    return [modified / (modified + unmodified) for modified, unmodified in observations]


def one_sided_mann_whitney(treat, control, minimum_replicates):
    if len(treat) < minimum_replicates or len(control) < minimum_replicates:
        return 1.0
    if treat and not isinstance(treat[0], (tuple, list, np.ndarray)):
        treat_rates = treat
        control_rates = control
    else:
        treat_rates = observation_rates(treat)
        control_rates = observation_rates(control)
    return float(
        mannwhitneyu(
            treat_rates,
            control_rates,
            alternative="greater",
            method="asymptotic",
        ).pvalue
    )


def _beta_binomial_log_likelihood(parameters, groups, separate_means):
    if separate_means:
        means = parameters[:2]
        rho = parameters[2]
    else:
        means = [parameters[0], parameters[0]]
        rho = parameters[1]

    concentration = (1.0 - rho) / rho
    total = 0.0
    for observations, mean in zip(groups, means):
        observations = np.asarray(observations, dtype=float)
        modified = observations[:, 0]
        unmodified = observations[:, 1]
        coverage = modified + unmodified
        alpha = mean * concentration
        beta = (1.0 - mean) * concentration
        total += float(
            np.sum(
                gammaln(coverage + 1.0)
                - gammaln(modified + 1.0)
                - gammaln(unmodified + 1.0)
                + betaln(modified + alpha, unmodified + beta)
                - betaln(alpha, beta)
            )
        )
    return total


def beta_binomial_lrt(treat, control, minimum_replicates):
    """One-sided treatment-greater beta-binomial likelihood-ratio test."""
    if len(treat) < minimum_replicates or len(control) < minimum_replicates:
        return 1.0

    def pooled_rate(observations):
        modified = sum(item[0] for item in observations)
        coverage = sum(sum(item) for item in observations)
        return modified / coverage if coverage else 0.0

    treat_rate = pooled_rate(treat)
    control_rate = pooled_rate(control)
    if treat_rate <= control_rate:
        return 1.0

    epsilon = 1e-6
    combined_rate = pooled_rate(treat + control)
    null = minimize(
        lambda values: -_beta_binomial_log_likelihood(values, (treat, control), False),
        x0=[np.clip(combined_rate, epsilon, 1.0 - epsilon), 0.05],
        bounds=[(epsilon, 1.0 - epsilon), (epsilon, 0.99)],
        method="L-BFGS-B",
    )
    alternative = minimize(
        lambda values: -_beta_binomial_log_likelihood(values, (treat, control), True),
        x0=[
            np.clip(treat_rate, epsilon, 1.0 - epsilon),
            np.clip(control_rate, epsilon, 1.0 - epsilon),
            0.05,
        ],
        bounds=[
            (epsilon, 1.0 - epsilon),
            (epsilon, 1.0 - epsilon),
            (epsilon, 0.99),
        ],
        method="L-BFGS-B",
    )
    if (
        not null.success
        or not alternative.success
        or not np.isfinite(null.fun)
        or not np.isfinite(alternative.fun)
    ):
        return 1.0
    statistic = max(0.0, 2.0 * (null.fun - alternative.fun))
    return float(min(1.0, 0.5 * chi2.sf(statistic, 1)))


def replicate_pvalue(treat, control, minimum_replicates, statistical_test):
    if statistical_test == "mannwhitney":
        return one_sided_mann_whitney(treat, control, minimum_replicates)
    if statistical_test == "beta_binomial":
        return beta_binomial_lrt(treat, control, minimum_replicates)
    raise ValueError(f"unknown replicate test: {statistical_test}")


def adjust_available(pvalues, available):
    adjusted = np.ones(len(pvalues), dtype=float)
    indexes = np.flatnonzero(available)
    if len(indexes):
        adjusted[indexes] = multipletests(
            np.asarray(pvalues)[indexes], method="fdr_bh"
        )[1]
    return adjusted


def confidence_status(candidate_eligible, replicate_available, replicate_pass):
    """Describe evidence strength without treating missing evidence as rejection."""
    if not candidate_eligible:
        return CONFIDENCE_STATUSES["excluded"]
    if not replicate_available:
        return CONFIDENCE_STATUSES["insufficient"]
    if replicate_pass:
        return CONFIDENCE_STATUSES["supported"]
    return CONFIDENCE_STATUSES["not_supported"]


def run_replicate_contrast(
    audit_path,
    calls_dir,
    samples_tsv,
    biotype,
    out_path,
    all_out_path=None,
    minimum_replicates=3,
    minimum_delta=0.1,
    maximum_padj=0.05,
    minimum_control_coverage=5,
    bam_dir=None,
    minimum_base_quality=20,
    minimum_mapping_quality=10,
    minimum_sample_coverage=1,
    statistical_test="beta_binomial",
    prefilter_effect=False,
    require_no_treat_significance=True,
    cell_line=None,
):
    audit = pd.read_csv(
        audit_path, sep="\t", dtype={"chrom": "string"}, low_memory=False
    )
    if audit.empty:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        audit.to_csv(out_path, sep="\t", index=False)
        if all_out_path:
            Path(all_out_path).parent.mkdir(parents=True, exist_ok=True)
            audit.to_csv(all_out_path, sep="\t", index=False)
        return audit, audit

    audit["start"] = pd.to_numeric(audit["start"], errors="raise").astype(int)
    coverage_mask = (
        (audit["coverage"] >= minimum_control_coverage)
        & (audit["pb_ctrl_coverage"] >= minimum_control_coverage)
        & (audit["no_treat_coverage"] >= minimum_control_coverage)
    )
    if prefilter_effect:
        candidate_mask = (
            coverage_mask
            & (audit["delta_pb_ctrl"] >= minimum_delta)
            & (audit["delta_no_treat"] >= minimum_delta)
            & (audit["mod_rate"] >= minimum_delta)
        )
    else:
        candidate_mask = coverage_mask
    audit["candidate_eligible"] = candidate_mask
    candidates = audit[candidate_mask].copy()
    candidates["site_key"] = list(
        zip(
            candidates["chrom"].astype(str),
            candidates["start"],
            candidates["strand"].astype(str),
        )
    )

    if candidates.empty:
        audit["replicate_pass"] = False
        audit["confidence_status"] = CONFIDENCE_STATUSES["excluded"]
        passing = audit.iloc[0:0].copy()
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        passing.to_csv(out_path, sep="\t", index=False)
        if all_out_path:
            Path(all_out_path).parent.mkdir(parents=True, exist_ok=True)
            audit.to_csv(all_out_path, sep="\t", index=False)
        return candidates, passing

    conditions = load_samples(samples_tsv, cell_line=cell_line)
    if bam_dir:
        rates = load_replicate_rates_from_bams(
            bam_dir,
            biotype,
            conditions,
            set(candidates["site_key"]),
            minimum_base_quality,
            minimum_mapping_quality,
            minimum_sample_coverage,
        )
    else:
        rates = load_replicate_rates(
            calls_dir, biotype, conditions, set(candidates["site_key"])
        )

    rows = []
    for row in candidates.itertuples(index=False):
        key = row.site_key
        groups = rates.get(key, {})
        treat = groups.get("treat", [])
        pb_ctrl = groups.get("pb_ctrl", [])
        no_treat = groups.get("no_treat", [])
        pb_available = (
            len(treat) >= minimum_replicates and len(pb_ctrl) >= minimum_replicates
        )
        untreated_available = (
            len(treat) >= minimum_replicates and len(no_treat) >= minimum_replicates
        )
        treat_median = float(np.median(observation_rates(treat))) if treat else np.nan
        pb_median = float(np.median(observation_rates(pb_ctrl))) if pb_ctrl else np.nan
        untreated_median = (
            float(np.median(observation_rates(no_treat))) if no_treat else np.nan
        )
        rows.append(
            {
                "treat_replicates": len(treat),
                "pb_ctrl_replicates": len(pb_ctrl),
                "no_treat_replicates": len(no_treat),
                "treat_median_rate": treat_median,
                "pb_ctrl_median_rate": pb_median,
                "no_treat_median_rate": untreated_median,
                "replicate_delta_pb_ctrl": treat_median - pb_median,
                "replicate_delta_no_treat": treat_median - untreated_median,
                "replicate_pb_pvalue": replicate_pvalue(
                    treat, pb_ctrl, minimum_replicates, statistical_test
                ),
                "replicate_no_treat_pvalue": replicate_pvalue(
                    treat, no_treat, minimum_replicates, statistical_test
                ),
                "replicate_test": statistical_test,
                "pb_available": pb_available,
                "no_treat_available": untreated_available,
            }
        )

    statistics = pd.DataFrame(rows, index=candidates.index)
    candidates = pd.concat([candidates, statistics], axis=1)
    candidates["replicate_pb_padj"] = adjust_available(
        candidates["replicate_pb_pvalue"], candidates["pb_available"]
    )
    candidates["replicate_no_treat_padj"] = adjust_available(
        candidates["replicate_no_treat_pvalue"],
        candidates["no_treat_available"],
    )
    significance_mask = candidates["replicate_pb_padj"] <= maximum_padj
    if require_no_treat_significance:
        significance_mask &= candidates["replicate_no_treat_padj"] <= maximum_padj
    candidates["replicate_pass"] = (
        candidates["pb_available"]
        & candidates["no_treat_available"]
        & (candidates["replicate_delta_pb_ctrl"] >= minimum_delta)
        & (candidates["replicate_delta_no_treat"] >= minimum_delta)
        & significance_mask
    )
    candidates["confidence_status"] = [
        confidence_status(True, pb and untreated, passed)
        for pb, untreated, passed in zip(
            candidates["pb_available"],
            candidates["no_treat_available"],
            candidates["replicate_pass"],
        )
    ]

    candidates["pooled_pvalue"] = candidates["pvalue"]
    candidates["pooled_padj"] = candidates["padj"]
    candidates["pvalue"] = candidates["replicate_pb_pvalue"]
    candidates["padj"] = candidates["replicate_pb_padj"]
    candidates = candidates.drop(columns=["site_key"])
    passing = candidates[candidates["replicate_pass"]].copy()

    audit_output = audit.merge(
        candidates.drop(
            columns=[column for column in audit.columns if column in candidates.columns]
        ),
        left_index=True,
        right_index=True,
        how="left",
    )
    audit_output["replicate_pass"] = audit_output["replicate_pass"].fillna(False)
    audit_output["confidence_status"] = audit_output.apply(
        lambda row: confidence_status(
            bool(row["candidate_eligible"]),
            bool(row.get("pb_available", False))
            and bool(row.get("no_treat_available", False)),
            bool(row["replicate_pass"]),
        ),
        axis=1,
    )

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    passing.to_csv(out_path, sep="\t", index=False)
    if all_out_path:
        Path(all_out_path).parent.mkdir(parents=True, exist_ok=True)
        audit_output.to_csv(all_out_path, sep="\t", index=False)
    return candidates, passing


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--calls-dir", required=True)
    parser.add_argument("--samples-tsv", required=True)
    parser.add_argument("--biotype", required=True)
    parser.add_argument("--cell-line")
    parser.add_argument("--out", required=True)
    parser.add_argument("--all-out")
    parser.add_argument("--minimum-replicates", type=int, default=3)
    parser.add_argument("--minimum-delta", type=float, default=0.1)
    parser.add_argument("--maximum-padj", type=float, default=0.05)
    parser.add_argument("--minimum-control-coverage", type=float, default=5)
    parser.add_argument(
        "--bam-dir",
        help="biotype BAM root; when set, derive replicate rates directly from BAMs",
    )
    parser.add_argument("--minimum-base-quality", type=int, default=20)
    parser.add_argument("--minimum-mapping-quality", type=int, default=10)
    parser.add_argument("--minimum-sample-coverage", type=float, default=1)
    parser.add_argument(
        "--statistical-test",
        choices=["beta_binomial", "mannwhitney"],
        default="beta_binomial",
    )
    parser.add_argument(
        "--prefilter-effect",
        action="store_true",
        help="legacy mode: select on pooled effect before replicate testing",
    )
    parser.add_argument(
        "--allow-nonsignificant-no-treat",
        action="store_true",
        help="do not require replicate significance against no-treatment",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tested, passing = run_replicate_contrast(
        audit_path=args.audit,
        calls_dir=args.calls_dir,
        samples_tsv=args.samples_tsv,
        biotype=args.biotype,
        out_path=args.out,
        all_out_path=args.all_out,
        minimum_replicates=args.minimum_replicates,
        minimum_delta=args.minimum_delta,
        maximum_padj=args.maximum_padj,
        minimum_control_coverage=args.minimum_control_coverage,
        bam_dir=args.bam_dir,
        minimum_base_quality=args.minimum_base_quality,
        minimum_mapping_quality=args.minimum_mapping_quality,
        minimum_sample_coverage=args.minimum_sample_coverage,
        statistical_test=args.statistical_test,
        prefilter_effect=args.prefilter_effect,
        require_no_treat_significance=not args.allow_nonsignificant_no_treat,
        cell_line=args.cell_line,
    )
    print(f"Replicate-tested sites: {len(tested):,}")
    print(f"Replicate-passing sites: {len(passing):,}")


if __name__ == "__main__":
    main()
