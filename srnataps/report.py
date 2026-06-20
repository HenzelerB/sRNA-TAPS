# -*- coding: utf-8 -*-
# Part of the sRNA-TAPS package (srnataps.report)
"""
srnataps.report
Generate an interactive HTML report from sRNA-TAPS pipeline outputs.

Uses Plotly for interactive figures and Jinja2 for HTML templating.
Output is a single self-contained HTML file — no server needed.

Usage:
    python -m srnataps.report --outdir /path/to/project
    srnataps module report --configfile config.yaml
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
from jinja2 import Environment, BaseLoader

from srnataps.utils import (
    detect_condition as get_condition,
    detect_cell_line as get_cell_line,
)

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s",
                    datefmt="%H:%M:%S")

# ── Colour schemes (matching R scripts) ──────────────────────────────────────
CONDITION_COLOURS = {
    "treat":    "#FFD680",  # Aurora gold  — TET + PB
    "pb_ctrl":  "#5BAFD0",  # Aurora blue  — PB only
    "no_treat": "#4DDEB8",  # Aurora mint  — untreated
}
CONDITION_LABELS = {
    "no_treat": "Untreated",
    "pb_ctrl":  "PB only",
    "treat":    "TET + PB",
}
BIOTYPE_COLOURS = {
    "miRNA":  "#FFD680",  # Aurora spectrum (R colorRampPalette, matches 00_setup.R)
    "tRNA":   "#D0CA96",
    "rRNA":   "#A1BFAD",
    "snoRNA": "#72B4C4",
    "snRNA":  "#59B5CC",
    "piRNA":  "#55C3C5",
    "lncRNA": "#51D0BE",
    "other":  "#4DDEB8",
}
TOOL_COLOURS = {
    "sRNA-TAPS": "#1C4062",  # brand navy
    "rastair":   "#CC79A7",  # Okabe-Ito reddish-purple
    "astair":    "#009E73",  # Okabe-Ito bluish-green
    "bismark":   "#E69F00",  # Okabe-Ito orange
}


# ══════════════════════════════════════════════════════════════════════════════
# Section builders — each returns a Plotly figure as JSON string
# ══════════════════════════════════════════════════════════════════════════════

def build_biotype_composition(outdir: str) -> str:
    """Stacked bar chart of biotype composition across all samples."""
    bio_file = Path(outdir) / "05.biotype_bams" / "biotype_composition_all_samples.tsv"
    if not bio_file.exists():
        return None

    df = pd.read_csv(bio_file, sep="\t")
    df["condition_label"] = df["condition"].map(CONDITION_LABELS).fillna(df["condition"])
    df["colour"]          = df["biotype"].map(BIOTYPE_COLOURS).fillna("#999999")

    # Sort samples
    df = df.sort_values(["cell_line", "condition", "sample"])

    fig = go.Figure()
    for biotype in list(BIOTYPE_COLOURS.keys()):
        sub = df[df["biotype"] == biotype]
        if len(sub) == 0:
            continue
        fig.add_trace(go.Bar(
            x    = sub["sample"],
            y    = sub["percent"],
            name = biotype,
            marker_color = BIOTYPE_COLOURS[biotype],
            hovertemplate = "<b>%{x}</b><br>" + biotype + ": %{y:.1f}%<extra></extra>",
        ))

    fig.update_layout(
        barmode     = "stack",
        title       = "RNA Biotype Composition",
        xaxis_title = "",
        yaxis_title = "% of mapped reads",
        legend_title= "Biotype",
        height      = 450,
        xaxis       = dict(tickangle=-45, tickfont=dict(size=9)),
        margin      = dict(b=120),
        plot_bgcolor= "white",
        paper_bgcolor="white",
    )
    fig.update_yaxes(range=[0, 100])
    return pio.to_json(fig)


def build_modrate_distribution(outdir: str, min_cov: int = 5) -> str:
    """Box plots of modification rate distribution per biotype per condition."""
    calls_dir = Path(outdir) / "07.taps_calls"
    if not calls_dir.exists():
        return None

    dfs = []
    for biotype in BIOTYPE_COLOURS:
        bt_dir = calls_dir / biotype
        if not bt_dir.exists():
            continue
        for f in bt_dir.glob("*_taps.tsv"):
            sample = f.name.replace(f"_{biotype}_taps.tsv", "")
            try:
                df = pd.read_csv(f, sep="\t")
                df = df[(df["coverage"] >= min_cov) & (df["snp_flag"] == "PASS")]
                if len(df) == 0:
                    continue
                df["biotype"]   = biotype
                df["sample"]    = sample
                df["condition"] = get_condition(sample)
                df["cell_line"] = get_cell_line(sample)
                dfs.append(df[["mod_rate", "biotype", "condition", "cell_line", "sample"]])
            except Exception:
                continue

    if not dfs:
        return None

    all_data = pd.concat(dfs, ignore_index=True)
    all_data = all_data[all_data["condition"].isin(["no_treat", "pb_ctrl", "treat"])]

    fig = go.Figure()
    for condition, label in CONDITION_LABELS.items():
        if condition not in ["no_treat", "pb_ctrl", "treat"]:
            continue
        sub = all_data[all_data["condition"] == condition]
        if len(sub) == 0:
            continue
        for biotype in BIOTYPE_COLOURS:
            bt_sub = sub[sub["biotype"] == biotype]
            if len(bt_sub) < 3:
                continue
            fig.add_trace(go.Box(
                y           = bt_sub["mod_rate"],
                name        = f"{label} | {biotype}",
                boxpoints   = False,
                marker_color= CONDITION_COLOURS[condition],
                legendgroup = condition,
                showlegend  = biotype == list(BIOTYPE_COLOURS.keys())[0],
                hovertemplate = f"<b>{label} | {biotype}</b><br>"
                                "Median: %{median:.3f}<extra></extra>",
            ))

    fig.update_layout(
        title       = "TAPS Modification Rate Distribution",
        yaxis_title = "Modification rate",
        xaxis_title = "",
        height      = 500,
        boxmode     = "group",
        plot_bgcolor= "white",
        paper_bgcolor="white",
        yaxis       = dict(range=[0, 1], tickformat=".0%"),
    )
    return pio.to_json(fig)


def build_top_species(outdir: str, min_cov: int = 5, top_n: int = 30) -> str:
    """Dot plot of top modified RNA species in TET+PB condition."""
    calls_dir = Path(outdir) / "07.taps_calls"
    if not calls_dir.exists():
        return None

    dfs = []
    for biotype in ["miRNA", "tRNA", "rRNA", "snoRNA"]:
        bt_dir = calls_dir / biotype
        if not bt_dir.exists():
            continue
        for f in bt_dir.glob("*treat*_taps.tsv"):
            sample = f.name.replace(f"_{biotype}_taps.tsv", "")
            if "no-treat" in sample or "pb_Ctrl" in sample:
                continue
            try:
                df = pd.read_csv(f, sep="\t")
                df = df[(df["coverage"] >= min_cov) & (df["snp_flag"] == "PASS")]
                df["biotype"]   = biotype
                df["sample"]    = sample
                df["cell_line"] = get_cell_line(sample)
                dfs.append(df)
            except Exception:
                continue

    if not dfs:
        return None

    all_data = pd.concat(dfs, ignore_index=True)
    if "gene_id" not in all_data.columns:
        return None

    summary = (
        all_data[all_data["gene_id"] != "."]
        .groupby(["biotype", "gene_id", "cell_line"])
        .agg(median_mod=("mod_rate", "median"),
             mean_cov=("coverage", "mean"),
             n_sites=("mod_rate", "count"))
        .reset_index()
        .query("n_sites >= 2")
    )

    top = (summary.sort_values("median_mod", ascending=False)
                  .groupby(["biotype", "cell_line"])
                  .head(top_n))

    fig = px.scatter(
        top,
        x         = "median_mod",
        y         = "gene_id",
        color     = "biotype",
        size      = "n_sites",
        facet_col = "cell_line",
        facet_row = "biotype",
        color_discrete_map = BIOTYPE_COLOURS,
        labels    = {"median_mod": "Median mod rate",
                     "gene_id": "", "n_sites": "N sites"},
        title     = f"Top {top_n} Modified RNA Species (TET+PB)",
        hover_data= ["mean_cov", "n_sites"],
    )
    fig.update_layout(
        height        = max(600, len(top) * 12),
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        showlegend    = False,
    )
    fig.update_xaxes(tickformat=".0%", range=[0, 1])
    return pio.to_json(fig)


def build_condition_scatter(outdir: str, min_cov: int = 5) -> str:
    """Scatter: PB-only vs TET+PB median mod_rate per gene."""
    calls_dir = Path(outdir) / "07.taps_calls"
    if not calls_dir.exists():
        return None

    dfs = {"pb_ctrl": [], "treat": []}
    for biotype in BIOTYPE_COLOURS:
        bt_dir = calls_dir / biotype
        if not bt_dir.exists():
            continue
        for f in bt_dir.glob("*_taps.tsv"):
            sample = f.name.replace(f"_{biotype}_taps.tsv", "")
            cond   = get_condition(sample)
            if cond not in dfs:
                continue
            try:
                df = pd.read_csv(f, sep="\t")
                df = df[(df["coverage"] >= min_cov) & (df["snp_flag"] == "PASS")]
                df["biotype"]   = biotype
                df["cell_line"] = get_cell_line(sample)
                dfs[cond].append(df)
            except Exception:
                continue

    if not dfs["pb_ctrl"] or not dfs["treat"]:
        return None

    def summarise(lst):
        return (pd.concat(lst)
                .query("gene_id != '.'")
                .groupby(["biotype", "gene_id", "cell_line"])["mod_rate"]
                .median().reset_index())

    pb    = summarise(dfs["pb_ctrl"]).rename(columns={"mod_rate": "pb_ctrl"})
    treat = summarise(dfs["treat"]).rename(columns={"mod_rate": "treat"})
    merged = pb.merge(treat, on=["biotype", "gene_id", "cell_line"])

    fig = px.scatter(
        merged,
        x             = "pb_ctrl",
        y             = "treat",
        color         = "biotype",
        facet_col     = "cell_line",
        hover_name    = "gene_id",
        color_discrete_map = BIOTYPE_COLOURS,
        labels        = {"pb_ctrl": "PB-only mod rate",
                         "treat":   "TET+PB mod rate"},
        title         = "Condition Comparison: PB-only vs TET+PB",
        opacity       = 0.7,
    )
    # Add y=x diagonal
    fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1,
                  line=dict(dash="dash", color="grey", width=1))
    fig.update_layout(
        height        = 500,
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
    )
    fig.update_xaxes(tickformat=".0%", range=[0, 1])
    fig.update_yaxes(tickformat=".0%", range=[0, 1])
    return pio.to_json(fig)


def build_concordance_heatmap(outdir: str) -> str:
    """Heatmap of Jaccard index between custom pipeline and benchmark tools."""
    conc_file = Path(outdir) / "09.compare" / "concordance_summary.tsv"
    if not conc_file.exists():
        return None

    df = pd.read_csv(conc_file, sep="\t")
    df = df[df["condition"] == "treat"]

    pivot = df.pivot_table(index="biotype", columns="tool",
                           values="jaccard", aggfunc="mean")

    fig = go.Figure(go.Heatmap(
        z           = pivot.values,
        x           = [t.capitalize() for t in pivot.columns],
        y           = pivot.index.tolist(),
        colorscale  = "YlOrRd",
        zmin        = 0, zmax = 1,
        text        = np.round(pivot.values, 2),
        texttemplate= "%{text}",
        hovertemplate = "Biotype: %{y}<br>Tool: %{x}<br>Jaccard: %{z:.3f}<extra></extra>",
    ))
    fig.update_layout(
        title         = "Concordance: Custom Pipeline vs Benchmark Tools (TET+PB)",
        xaxis_title   = "Benchmark tool",
        yaxis_title   = "RNA biotype",
        height        = 450,
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
    )
    return pio.to_json(fig)


def build_correlation_plot(outdir: str) -> str:
    """Line plot of Pearson correlation per biotype per tool."""
    corr_file = Path(outdir) / "09.compare" / "correlation_summary.tsv"
    if not corr_file.exists():
        return None

    df = pd.read_csv(corr_file, sep="\t")
    df = df[(df["condition"] == "treat") & df["pearson_r"].notna()]

    fig = go.Figure()
    for tool, colour in TOOL_COLOURS.items():
        sub = df[df["tool"] == tool]
        if len(sub) == 0:
            continue
        fig.add_trace(go.Scatter(
            x    = sub["biotype"],
            y    = sub["pearson_r"],
            name = tool.capitalize(),
            mode = "lines+markers",
            line = dict(color=colour, width=2),
            marker = dict(size=8),
            hovertemplate = f"<b>{tool}</b><br>%{{x}}<br>Pearson r = %{{y:.3f}}<extra></extra>",
        ))

    fig.add_hline(y=0, line_dash="dash", line_color="grey", line_width=1)
    fig.update_layout(
        title         = "Pearson Correlation at Shared Sites (TET+PB)",
        xaxis_title   = "RNA biotype",
        yaxis_title   = "Pearson r",
        height        = 400,
        yaxis         = dict(range=[-1, 1]),
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
    )
    return pio.to_json(fig)


# ══════════════════════════════════════════════════════════════════════════════
# HTML template
# ══════════════════════════════════════════════════════════════════════════════

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>sRNA-TAPS Report — {{ project_name }}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
         background: #f8f9fa; color: #212529; }
  .header { background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
            color: white; padding: 32px 40px; }
  .header h1 { font-size: 1.8rem; font-weight: 700; }
  .header p  { margin-top: 6px; opacity: 0.85; font-size: 0.95rem; }
  .nav { background: white; border-bottom: 1px solid #dee2e6;
         padding: 0 40px; display: flex; gap: 0; position: sticky; top: 0; z-index: 100; }
  .nav a { display: inline-block; padding: 14px 20px; text-decoration: none;
           color: #495057; font-size: 0.88rem; font-weight: 500; border-bottom: 3px solid transparent; }
  .nav a:hover, .nav a.active { color: #1a237e; border-bottom-color: #1a237e; }
  .content { max-width: 1400px; margin: 0 auto; padding: 32px 40px; }
  .section { margin-bottom: 48px; }
  .section h2 { font-size: 1.3rem; font-weight: 700; color: #1a237e;
                border-left: 4px solid #1a237e; padding-left: 12px; margin-bottom: 20px; }
  .section h3 { font-size: 1rem; font-weight: 600; color: #495057; margin: 24px 0 10px; }
  .card { background: white; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
          padding: 20px; margin-bottom: 20px; }
  .plot-container { width: 100%; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
                gap: 16px; margin-bottom: 24px; }
  .stat-card { background: white; border-radius: 8px; padding: 16px 20px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.08); border-top: 3px solid #1a237e; }
  .stat-card .value { font-size: 1.6rem; font-weight: 700; color: #1a237e; }
  .stat-card .label { font-size: 0.78rem; color: #6c757d; margin-top: 4px; }
  .note { background: #e8f4fd; border-left: 4px solid #2196F3; padding: 12px 16px;
          border-radius: 0 6px 6px 0; font-size: 0.85rem; color: #1565C0; margin: 12px 0; }
  .footer { text-align: center; padding: 32px; color: #6c757d; font-size: 0.82rem;
            border-top: 1px solid #dee2e6; margin-top: 48px; }
  @media (max-width: 768px) { .content { padding: 16px; } .nav { overflow-x: auto; } }
</style>
</head>
<body>

<div class="header">
  <h1>sRNA-TAPS Analysis Report</h1>
  <p>Project: <strong>{{ project_name }}</strong> &nbsp;|&nbsp;
     Generated: {{ timestamp }} &nbsp;|&nbsp;
     sRNA-TAPS v{{ version }}</p>
</div>

<nav class="nav">
  <a href="#summary">Summary</a>
  <a href="#qc">QC</a>
  <a href="#biotype">Biotypes</a>
  <a href="#modification">Modifications</a>
  <a href="#benchmarking">Benchmarking</a>
</nav>

<div class="content">

  <!-- Summary stats -->
  <div class="section" id="summary">
    <h2>Summary</h2>
    <div class="stats-grid">
      <div class="stat-card">
        <div class="value">{{ n_samples }}</div>
        <div class="label">Total samples</div>
      </div>
      <div class="stat-card">
        <div class="value">{{ n_conditions }}</div>
        <div class="label">Conditions</div>
      </div>
      <div class="stat-card">
        <div class="value">{{ n_cell_lines }}</div>
        <div class="label">Cell lines</div>
      </div>
      <div class="stat-card">
        <div class="value">{{ n_sites_total }}</div>
        <div class="label">Total PASS sites</div>
      </div>
      <div class="stat-card">
        <div class="value">{{ n_sites_treat }}</div>
        <div class="label">Sites in TET+PB</div>
      </div>
    </div>
    <div class="note">
      TAPS chemistry: m5C and 5hmC → T (TET oxidation + pyridine borane reduction).
      Unmodified C stays as C. C→T in reads = modification signal (opposite of bisulfite).
      Three-layer SNP filtering applied before calling.
    </div>
  </div>

  <!-- Biotype composition -->
  <div class="section" id="biotype">
    <h2>RNA Biotype Composition</h2>
    <div class="card">
      {% if biotype_fig %}
      <div id="biotype-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Biotype composition data not available.</p>
      {% endif %}
    </div>
  </div>

  <!-- Modification rates -->
  <div class="section" id="modification">
    <h2>TAPS Modification Analysis</h2>

    <h3>Modification rate distribution per biotype</h3>
    <div class="card">
      {% if modrate_fig %}
      <div id="modrate-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Modification rate data not available.</p>
      {% endif %}
    </div>

    <h3>Top modified RNA species (TET+PB)</h3>
    <div class="card">
      {% if topspecies_fig %}
      <div id="topspecies-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Top species data not available.</p>
      {% endif %}
    </div>

    <h3>Condition comparison: PB-only vs TET+PB</h3>
    <div class="card">
      <div class="note">
        Sites above the diagonal are enriched by TET oxidation — these are genuine
        m5C or 5hmC sites. Sites on the diagonal are chemical background from
        pyridine borane alone.
      </div>
      {% if scatter_fig %}
      <div id="scatter-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Condition comparison data not available.</p>
      {% endif %}
    </div>
  </div>

  <!-- Benchmarking -->
  <div class="section" id="benchmarking">
    <h2>Benchmarking Comparison</h2>
    <div class="note">
      All three tools use Bowtie1-aligned BAMs as input (aligner matched).
      Bismark chemistry is inverted (1 − rate) before comparison — negative
      correlation before inversion is expected and serves as a positive control.
    </div>

    <h3>Site-level concordance (Jaccard index)</h3>
    <div class="card">
      {% if concordance_fig %}
      <div id="concordance-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Benchmarking data not available. Run with --benchmark flag.</p>
      {% endif %}
    </div>

    <h3>Pearson correlation at shared sites</h3>
    <div class="card">
      {% if correlation_fig %}
      <div id="correlation-plot" class="plot-container"></div>
      {% else %}
      <p style="color:#6c757d">Correlation data not available.</p>
      {% endif %}
    </div>
  </div>

</div>

<div class="footer">
  Generated by sRNA-TAPS v{{ version }} on {{ timestamp }}<br>
  <a href="https://github.com/HenzelerB/sRNA-TAPS" style="color:#1a237e">
    github.com/HenzelerB/sRNA-TAPS
  </a>
</div>

<script>
{% if biotype_fig %}
Plotly.newPlot('biotype-plot', JSON.parse('{{ biotype_fig | tojson }}').data,
               JSON.parse('{{ biotype_fig | tojson }}').layout, {responsive: true});
{% endif %}
{% if modrate_fig %}
Plotly.newPlot('modrate-plot', JSON.parse('{{ modrate_fig | tojson }}').data,
               JSON.parse('{{ modrate_fig | tojson }}').layout, {responsive: true});
{% endif %}
{% if topspecies_fig %}
Plotly.newPlot('topspecies-plot', JSON.parse('{{ topspecies_fig | tojson }}').data,
               JSON.parse('{{ topspecies_fig | tojson }}').layout, {responsive: true});
{% endif %}
{% if scatter_fig %}
Plotly.newPlot('scatter-plot', JSON.parse('{{ scatter_fig | tojson }}').data,
               JSON.parse('{{ scatter_fig | tojson }}').layout, {responsive: true});
{% endif %}
{% if concordance_fig %}
Plotly.newPlot('concordance-plot', JSON.parse('{{ concordance_fig | tojson }}').data,
               JSON.parse('{{ concordance_fig | tojson }}').layout, {responsive: true});
{% endif %}
{% if correlation_fig %}
Plotly.newPlot('correlation-plot', JSON.parse('{{ correlation_fig | tojson }}').data,
               JSON.parse('{{ correlation_fig | tojson }}').layout, {responsive: true});
{% endif %}
</script>
</body>
</html>
"""


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir",   required=True, help="Project output directory")
    p.add_argument("--out",      default=None,  help="Output HTML path (default: outdir/report/srnataps_report.html)")
    p.add_argument("--min-cov",  type=int, default=5)
    p.add_argument("--project",  default=None,  help="Project name for report header")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = args.outdir

    report_dir = Path(outdir) / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    out_html = args.out or str(report_dir / "srnataps_report.html")
    project  = args.project or Path(outdir).name

    log.info("Building interactive HTML report...")
    log.info("  Project : %s", project)
    log.info("  Output  : %s", out_html)

    # ── Build all figures ─────────────────────────────────────────────────────
    log.info("Building figures...")

    biotype_fig    = build_biotype_composition(outdir)
    log.info("  Biotype composition : %s", "OK" if biotype_fig else "SKIP")

    modrate_fig    = build_modrate_distribution(outdir, args.min_cov)
    log.info("  Mod rate dist       : %s", "OK" if modrate_fig else "SKIP")

    topspecies_fig = build_top_species(outdir, args.min_cov)
    log.info("  Top species         : %s", "OK" if topspecies_fig else "SKIP")

    scatter_fig    = build_condition_scatter(outdir, args.min_cov)
    log.info("  Condition scatter   : %s", "OK" if scatter_fig else "SKIP")

    concordance_fig = build_concordance_heatmap(outdir)
    log.info("  Concordance heatmap : %s", "OK" if concordance_fig else "SKIP")

    correlation_fig = build_correlation_plot(outdir)
    log.info("  Correlation plot    : %s", "OK" if correlation_fig else "SKIP")

    # ── Count summary stats ───────────────────────────────────────────────────
    calls_dir = Path(outdir) / "07.taps_calls"
    n_sites_total = n_sites_treat = "N/A"
    if calls_dir.exists():
        try:
            total = sum(
                len(pd.read_csv(f, sep="\t").query("snp_flag == 'PASS'"))
                for f in calls_dir.rglob("*_taps.tsv")
            )
            treat = sum(
                len(pd.read_csv(f, sep="\t").query("snp_flag == 'PASS'"))
                for f in calls_dir.rglob("*treat*_taps.tsv")
                if "no-treat" not in str(f) and "pb_Ctrl" not in str(f)
            )
            n_sites_total = f"{total:,}"
            n_sites_treat = f"{treat:,}"
        except Exception:
            pass

    # ── Render HTML ───────────────────────────────────────────────────────────
    from srnataps import __version__

    env      = Environment(loader=BaseLoader())
    template = env.from_string(HTML_TEMPLATE)

    # Sample / condition / cell-line counts from the sample sheet (authoritative,
    # matches 00_setup.R). Falls back to "N/A" if the sheet is absent.
    n_samples = n_conditions = n_cell_lines = "N/A"
    samples_tsv = Path(outdir) / "samples.tsv"
    if samples_tsv.exists():
        try:
            _sdf = pd.read_csv(samples_tsv, sep="\t")
            n_samples    = str(len(_sdf))
            n_conditions = str(_sdf["condition"].nunique())
            n_cell_lines = str(_sdf["cell_line"].nunique())
        except Exception:
            pass

    html = template.render(
        project_name    = project,
        timestamp       = datetime.now().strftime("%Y-%m-%d %H:%M"),
        version         = __version__,
        n_samples       = n_samples,
        n_conditions    = n_conditions,
        n_cell_lines    = n_cell_lines,
        n_sites_total   = n_sites_total,
        n_sites_treat   = n_sites_treat,
        biotype_fig     = biotype_fig,
        modrate_fig     = modrate_fig,
        topspecies_fig  = topspecies_fig,
        scatter_fig     = scatter_fig,
        concordance_fig = concordance_fig,
        correlation_fig = correlation_fig,
    )

    with open(out_html, "w") as f:
        f.write(html)

    log.info("Report written: %s", out_html)


if __name__ == "__main__":
    main()
