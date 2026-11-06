"""Build a compact presentation bundle from completed sRNA-TAPS outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BIOTYPE_ORDER = ["miRNA", "tRNA", "rRNA", "snoRNA", "snRNA", "lncRNA", "other", "piRNA"]
BIOTYPE_COLORS = {
    "miRNA": "#FFD680",
    "tRNA": "#D0CA96",
    "rRNA": "#A1BFAD",
    "snoRNA": "#72B4C4",
    "snRNA": "#59B5CC",
    "lncRNA": "#51D0BE",
    "other": "#4DDEB8",
    "piRNA": "#A8CADA",
}


def _read_tsv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype={"chrom": str}, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def pooled_signal_summary(outdir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted((outdir / "07b.pooled_calls").glob("*/treat_*_pooled_*_taps.tsv")):
        biotype = path.parent.name
        name = path.name
        cell_line = name.removeprefix("treat_").split("_pooled_")[0]
        df = _read_tsv(path)
        if df.empty:
            rows.append(
                {
                    "cell_line": cell_line,
                    "biotype": biotype,
                    "sites": 0,
                    "rate_ge_0.01": 0,
                    "rate_ge_0.05": 0,
                    "rate_ge_0.10": 0,
                    "median_rate": 0.0,
                    "median_coverage": 0.0,
                }
            )
            continue
        rows.append(
            {
                "cell_line": cell_line,
                "biotype": biotype,
                "sites": len(df),
                "rate_ge_0.01": int((df["mod_rate"] >= 0.01).sum()),
                "rate_ge_0.05": int((df["mod_rate"] >= 0.05).sum()),
                "rate_ge_0.10": int((df["mod_rate"] >= 0.10).sum()),
                "median_rate": float(df["mod_rate"].median()),
                "median_coverage": float(df["coverage"].median()) if "coverage" in df else 0.0,
            }
        )
    return pd.DataFrame(rows)


def evidence_summary(outdir: Path) -> pd.DataFrame:
    path = outdir / "report_3letter" / "tables" / "evidence_tier_summary.tsv"
    if not path.exists():
        path = outdir / "report" / "tables" / "evidence_tier_summary.tsv"
    return _read_tsv(path) if path.exists() else pd.DataFrame()


def shared_gene_summary(outdir: Path) -> pd.DataFrame:
    path = outdir / "07x.no_treat_only_contrast" / "shared_HEK_Caco2_gene_summary.tsv"
    return _read_tsv(path) if path.exists() else pd.DataFrame()


def _plot_pooled(ax, pooled: pd.DataFrame) -> None:
    if pooled.empty:
        ax.text(0.5, 0.5, "No pooled call table found", ha="center", va="center")
        ax.axis("off")
        return
    plot = pooled[pooled["biotype"].isin(["miRNA", "tRNA", "rRNA", "snoRNA", "snRNA"])].copy()
    plot["biotype"] = pd.Categorical(plot["biotype"], BIOTYPE_ORDER, ordered=True)
    plot = plot.sort_values(["biotype", "cell_line"])
    pivot = plot.pivot_table(
        index="biotype", columns="cell_line", values="rate_ge_0.10", aggfunc="sum", fill_value=0
    )
    pivot.plot(kind="bar", ax=ax, color=["#06D6A0", "#FFD166"], width=0.78)
    ax.set_title("A. Treatment signal at pooled sites")
    ax.set_ylabel("Sites with mod_rate >= 10%")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Cell line", frameon=False)


def _plot_evidence(ax, evidence: pd.DataFrame) -> None:
    if evidence.empty:
        ax.text(0.5, 0.5, "No evidence tier table found", ha="center", va="center")
        ax.axis("off")
        return
    plot = evidence[evidence["biotype"].isin(["miRNA", "tRNA", "rRNA", "snoRNA", "snRNA"])].copy()
    if plot.empty:
        ax.text(0.5, 0.5, "No canonical sRNA evidence tiers found", ha="center", va="center")
        ax.axis("off")
        return
    tier_col = "evidence_tier" if "evidence_tier" in plot.columns else "tier"
    count_col = "sites" if "sites" in plot.columns else "n_sites"
    plot[tier_col] = pd.Categorical(
        plot[tier_col],
        ["pooled", "control_contrast", "replicate_discovery", "stringent"],
        ordered=True,
    )
    totals = plot.groupby(tier_col, observed=False)[count_col].sum().reindex(
        ["pooled", "control_contrast", "replicate_discovery", "stringent"], fill_value=0
    )
    ax.plot(range(len(totals)), totals.values, marker="o", color="#1C4062", linewidth=2)
    ax.fill_between(range(len(totals)), totals.values, color="#A8CADA", alpha=0.35)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(range(len(totals)))
    ax.set_xticklabels(["pooled", "dual ctrl", "replicate", "stringent"], rotation=25, ha="right")
    ax.set_ylabel("Total sites")
    ax.set_title("B. Current dual-control tiers")


def _plot_shared_genes(ax, genes: pd.DataFrame) -> None:
    if genes.empty:
        ax.text(0.5, 0.5, "No shared HEK/Caco2 no_treat-primary table found", ha="center", va="center")
        ax.axis("off")
        return
    name_col = "gene_name" if "gene_name" in genes.columns else "gene_id"
    count_col = "n_sites" if "n_sites" in genes.columns else "sites"
    plot = genes.copy()
    if "gene_biotype" in plot:
        plot = plot[plot["gene_biotype"].isin(["miRNA", "snoRNA", "rRNA", "Mt_tRNA", "snRNA", "tRNA"])]
    plot = plot.sort_values(count_col, ascending=False).head(12).iloc[::-1]
    colors = [BIOTYPE_COLORS.get(bt, "#7CBBD4") for bt in plot.get("gene_biotype", "")]
    ax.barh(plot[name_col].astype(str), plot[count_col], color=colors)
    ax.set_xlabel("Shared sites")
    ax.set_title("C. Recurrent no_treat-primary candidates")


def _plot_message(ax) -> None:
    ax.axis("off")
    message = (
        "D. Interpretation for presentation\n\n"
        "Primary signal: TET+PB treatment is compared to no_treat.\n\n"
        "PB-only is retained as a process-control/background layer, not as a hard "
        "subtraction gate, because PB can still permit conversion in this dataset.\n\n"
        "The dual-control stringent tier is therefore conservative here; the "
        "no_treat-primary/shared sites are the biologically useful discovery set."
    )
    ax.text(
        0.02,
        0.98,
        message,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#F5F8FA", "edgecolor": "#A8CADA"},
        wrap=True,
    )


def make_figure(outdir: Path, bundle_dir: Path, pooled: pd.DataFrame, evidence: pd.DataFrame, genes: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    _plot_pooled(axes[0, 0], pooled)
    _plot_evidence(axes[0, 1], evidence)
    _plot_shared_genes(axes[1, 0], genes)
    _plot_message(axes[1, 1])
    fig.suptitle("sRNA-TAPS bio3 summary: treatment signal and evidence interpretation", fontsize=13)
    bundle_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(bundle_dir / f"bio3_final_summary.{ext}", dpi=300)
    plt.close(fig)


def write_markdown(bundle_dir: Path, pooled: pd.DataFrame, evidence: pd.DataFrame, genes: pd.DataFrame) -> None:
    lines = [
        "# sRNA-TAPS bio3 presentation summary",
        "",
        "## Main conclusion",
        "",
        "PB-only behaves as a process-control/background layer in this dataset. "
        "The presentation should emphasize treatment vs no_treat as the primary "
        "discovery contrast, with PB reported alongside rather than used as a hard subtraction gate.",
        "",
        "## Generated files",
        "",
        "- `figures/bio3_final_summary.pdf`",
        "- `figures/bio3_final_summary.png`",
        "- `figures/bio3_final_summary.svg`",
        "- `tables/pooled_treatment_signal_summary.tsv`",
        "- `tables/evidence_tier_summary.tsv`",
        "- `tables/shared_HEK_Caco2_gene_summary.tsv` if available",
        "",
    ]
    if not pooled.empty:
        total = int(pooled["rate_ge_0.10"].sum())
        lines += [
            "## Pooled treatment signal",
            "",
            f"Treatment pooled calls contain {total:,} sites with modification rate >= 10% across reported biotypes/cell lines.",
            "",
        ]
    if not evidence.empty:
        tier_col = "evidence_tier" if "evidence_tier" in evidence.columns else "tier"
        count_col = "sites" if "sites" in evidence.columns else "n_sites"
        tier_counts = evidence.groupby(tier_col, observed=False)[count_col].sum()
        lines += ["## Current evidence tiers", ""]
        for tier, count in tier_counts.items():
            lines.append(f"- {tier}: {int(count):,} sites")
        lines.append("")
    if not genes.empty:
        name_col = "gene_name" if "gene_name" in genes.columns else "gene_id"
        count_col = "n_sites" if "n_sites" in genes.columns else "sites"
        top = genes.sort_values(count_col, ascending=False).head(10)
        lines += ["## Top shared no_treat-primary genes", ""]
        for row in top.itertuples():
            lines.append(f"- {getattr(row, name_col)}: {int(getattr(row, count_col))} shared sites")
        lines.append("")
    (bundle_dir / "bio3_presentation_summary.md").write_text("\n".join(lines), encoding="utf-8")


def build_bundle(outdir: Path, output: Path) -> None:
    bundle = output if output.is_absolute() else outdir / output
    fig_dir = bundle / "figures"
    table_dir = bundle / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    pooled = pooled_signal_summary(outdir)
    evidence = evidence_summary(outdir)
    genes = shared_gene_summary(outdir)

    if not pooled.empty:
        _write_tsv(pooled, table_dir / "pooled_treatment_signal_summary.tsv")
    if not evidence.empty:
        _write_tsv(evidence, table_dir / "evidence_tier_summary.tsv")
    if not genes.empty:
        _write_tsv(genes, table_dir / "shared_HEK_Caco2_gene_summary.tsv")

    make_figure(outdir, fig_dir, pooled, evidence, genes)
    write_markdown(bundle, pooled, evidence, genes)
    print(f"Presentation bundle written: {bundle}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, help="Completed sRNA-TAPS run directory")
    parser.add_argument(
        "--output",
        default="presentation_bundle",
        help="Output directory, relative to --outdir unless absolute",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_bundle(Path(args.outdir).resolve(), Path(args.output))


if __name__ == "__main__":
    main()
