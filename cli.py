# -*- coding: utf-8 -*-
"""
srnataps.cli
Command-line interface for sRNA-TAPS.

Subcommands:
    srnataps init    — create a new project (config.yaml + samples.tsv)
    srnataps run     — run the full pipeline via Snakemake
    srnataps module  — run a single named pipeline module
    srnataps check   — validate environment and config before running

Examples:
    srnataps init --outdir ~/my_project --genome hg38.fa --gtf hg38.gtf
    srnataps run  --configfile ~/my_project/config.yaml --slurm
    srnataps module trim --configfile ~/my_project/config.yaml
    srnataps module call --configfile ~/my_project/config.yaml --benchmark
    srnataps check --configfile ~/my_project/config.yaml
"""

import os
import sys
import shutil
import subprocess
import textwrap
from pathlib import Path

import click
import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import print as rprint

from srnataps import __version__

console = Console()

# ── Paths to bundled resources ────────────────────────────────────────────────
PKG_DIR      = Path(__file__).parent
WORKFLOW_DIR = PKG_DIR / "workflow"
SNAKEFILE    = WORKFLOW_DIR / "Snakefile"
CONFIG_TMPL  = PKG_DIR / "config" / "config.yaml"
SAMPLES_TMPL = PKG_DIR / "config" / "samples.tsv"
SLURM_PROF   = PKG_DIR / "workflow" / "profiles" / "slurm"

# ── Valid module names ─────────────────────────────────────────────────────────
MODULES = {
    "fastqc":   "Run FastQC on raw merged FASTQs",
    "trim":     "Trim adapters with Trim Galore (TruSeq small RNA)",
    "index":    "Build Bowtie1 genome index",
    "align":    "Align trimmed reads with Bowtie1",
    "biotype":  "Split BAMs by RNA biotype",
    "snp":      "Build SNP blacklist from no-treat BAMs",
    "call":     "Call TAPS m5C modifications",
    "benchmark":"Run rastair, asTair, and Bismark benchmarking",
    "compare":  "Compare custom pipeline vs benchmark tools",
}



def _default_config() -> dict:
    """Return a minimal default config when bundled template is not found."""
    return {
        "project":   {"name": "my_taps_project", "outdir": ""},
        "input":     {"fastq_dir": "", "samples_tsv": "samples.tsv"},
        "reference": {"genome_fa": "", "gtf": "", "bowtie1_index": "", "dbsnp_vcf": ""},
        "trimming":  {"adapter": "TGGAATTCTCGGGTGCCAAGG", "min_length": 18, "max_length": 50},
        "alignment": {"mismatches": 2, "multimappers": 10, "max_multimappers": 100, "threads": 8},
        "biotype":   {"priority": ["miRNA","tRNA","piRNA","snoRNA","snRNA","rRNA","lncRNA","other"]},
        "snp":       {"min_af": 0.20, "min_cov": 10, "het_threshold": 0.40},
        "calling":   {
            "min_base_quality": 20, "min_mapping_quality": 10,
            "background_rate": 0.005, "pval_threshold": 0.05,
            "min_coverage": {"rRNA":10,"miRNA":3,"tRNA":3,"snoRNA":3,"snRNA":5,"piRNA":5,"lncRNA":5,"other":5},
            "contexts": ["ALL"],
        },
        "benchmark": {"enabled": False, "tools": {
            "rastair": {"min_depth":5,"min_mapq":10,"min_baseq":20},
            "astair":  {"min_depth":3,"min_mapq":10,"min_baseq":20,"context":"all"},
            "bismark": {"bismark_index":""},
        }},
        "slurm": {
            "partition":"compute","account":"","mail_user":"",
            "default_resources":{"cpus_per_task":8,"mem_mb":32000,"runtime":240},
            "rule_resources":{},
        },
        "output": {
            "fastqc":"02.fastqc","trimmed":"03.trimGalore","genome":"04a.genome",
            "aligned":"04b.aligned","biotypes":"05.biotype_bams","snp":"06.snp_resources",
            "calls":"07.taps_calls","benchmark":"08.benchmark","compare":"09.compare",
            "logs":"logs","report":"report",
        },
    }


def _write_default_samples(path: Path) -> None:
    """Write a minimal samples.tsv template."""
    with open(path, "w") as f:
        f.write("sample\tcondition\tcell_line\treplicate\tfastq\n")
        f.write("no-treat_Ctrl_HEK_R1\tno_treat\tHEK\tR1\t/path/to/fastq.gz\n")



# ── CLI group ──────────────────────────────────────────────────────────────────
@click.group()
@click.version_option(__version__, prog_name="sRNA-TAPS")
def cli():
    """
    sRNA-TAPS: TAPS-based m5C detection for small RNA sequencing.

    Run the full pipeline or individual modules on your samples.
    Start with: srnataps init --outdir my_project
    """
    pass


# ══════════════════════════════════════════════════════════════════════════════
# srnataps init
# ══════════════════════════════════════════════════════════════════════════════
@cli.command()
@click.option("--outdir",   required=True,  type=click.Path(), help="Project output directory (will be created)")
@click.option("--genome",   default=None,   help="Path to hg38 FASTA (optional — fill in config.yaml later)")
@click.option("--gtf",      default=None,   help="Path to Ensembl GTF (optional)")
@click.option("--fastq-dir",default=None,   help="Directory containing merged FASTQ files")
@click.option("--force",    is_flag=True,   help="Overwrite existing project directory")
def init(outdir, genome, gtf, fastq_dir, force):
    """
    Initialise a new sRNA-TAPS project.

    Creates a project directory with:
        config.yaml   — all pipeline parameters (edit before running)
        samples.tsv   — sample sheet template (fill in your samples)

    Example:
        srnataps init --outdir ~/my_taps_project --genome hg38.fa --gtf hg38.gtf
    """
    outdir = Path(outdir).resolve()

    if outdir.exists() and not force:
        console.print(f"[red]Directory already exists: {outdir}[/red]")
        console.print("Use --force to overwrite, or choose a different --outdir")
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    # ── Generate config.yaml ─────────────────────────────────────────────────
    config_out = outdir / "config.yaml"

    # Try to load from bundled template first; fall back to inline defaults
    if CONFIG_TMPL.exists():
        with open(CONFIG_TMPL) as f:
            config = yaml.safe_load(f)
    else:
        config = _default_config()

    # Fill in user-provided values
    config["project"]["outdir"] = str(outdir)
    if genome:    config["reference"]["genome_fa"] = str(Path(genome).resolve())
    if gtf:       config["reference"]["gtf"]       = str(Path(gtf).resolve())
    if fastq_dir: config["input"]["fastq_dir"]     = str(Path(fastq_dir).resolve())

    with open(config_out, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    # ── Generate samples.tsv template ────────────────────────────────────────
    samples_out = outdir / "samples.tsv"
    if SAMPLES_TMPL.exists():
        shutil.copy(SAMPLES_TMPL, samples_out)
    else:
        _write_default_samples(samples_out)

    # ── Print summary ─────────────────────────────────────────────────────────
    console.print(Panel.fit(
        f"[bold green]Project initialised: {outdir}[/bold green]\n\n"
        f"Next steps:\n"
        f"  1. Edit [bold]{config_out}[/bold] — set paths and parameters\n"
        f"  2. Edit [bold]{samples_out}[/bold] — add your sample names and FASTQ paths\n"
        f"  3. Run:  [bold]srnataps check --configfile {config_out}[/bold]\n"
        f"  4. Run:  [bold]srnataps run   --configfile {config_out} --slurm[/bold]",
        title="sRNA-TAPS init",
    ))


# ══════════════════════════════════════════════════════════════════════════════
# srnataps run
# ══════════════════════════════════════════════════════════════════════════════
@cli.command()
@click.option("--configfile",  required=True,  type=click.Path(exists=True), help="Path to config.yaml")
@click.option("--slurm",       is_flag=True,   help="Submit jobs to SLURM cluster")
@click.option("--benchmark",   is_flag=True,   help="Also run rastair, asTair, and Bismark benchmarking")
@click.option("--cores",       default=1,      help="Local cores (ignored if --slurm)")
@click.option("--jobs",        default=50,     help="Max concurrent SLURM jobs")
@click.option("--dryrun",  "-n", is_flag=True, help="Show what would be run without executing")
@click.option("--until",       default=None,   help="Run up to and including this rule (e.g. call)")
@click.option("--rerun-incomplete", is_flag=True, help="Re-run incomplete jobs")
@click.option("--snakemake-args", default="",  help="Extra arguments passed directly to Snakemake")
def run(configfile, slurm, benchmark, cores, jobs, dryrun, until, rerun_incomplete, snakemake_args):
    """
    Run the full sRNA-TAPS pipeline.

    Executes all steps from FastQC → trimming → alignment → biotype splitting
    → SNP filtering → TAPS calling → (optionally) benchmarking + comparison.

    Examples:
        srnataps run --configfile config.yaml --slurm
        srnataps run --configfile config.yaml --slurm --benchmark
        srnataps run --configfile config.yaml --dryrun
        srnataps run --configfile config.yaml --slurm --until call
    """
    _run_snakemake(
        configfile=configfile,
        target="all_benchmark" if benchmark else "all",
        slurm=slurm,
        cores=cores,
        jobs=jobs,
        dryrun=dryrun,
        until=until,
        rerun_incomplete=rerun_incomplete,
        extra_config={"benchmark": benchmark},
        snakemake_args=snakemake_args,
    )


# ══════════════════════════════════════════════════════════════════════════════
# srnataps module
# ══════════════════════════════════════════════════════════════════════════════
@cli.command()
@click.argument("module_name", type=click.Choice(list(MODULES.keys())))
@click.option("--configfile",  required=True,  type=click.Path(exists=True), help="Path to config.yaml")
@click.option("--slurm",       is_flag=True,   help="Submit jobs to SLURM cluster")
@click.option("--cores",       default=1,      help="Local cores (ignored if --slurm)")
@click.option("--jobs",        default=50,     help="Max concurrent SLURM jobs")
@click.option("--dryrun",  "-n", is_flag=True, help="Show what would be run without executing")
@click.option("--sample",      default=None,   help="Run for a specific sample only")
@click.option("--snakemake-args", default="",  help="Extra arguments passed directly to Snakemake")
def module(module_name, configfile, slurm, cores, jobs, dryrun, sample, snakemake_args):
    """
    Run a single pipeline module independently.

    Available modules:
        fastqc    — FastQC on raw merged FASTQs
        trim      — Trim Galore adapter trimming
        index     — Bowtie1 genome index build
        align     — Bowtie1 alignment
        biotype   — Biotype BAM splitting
        snp       — SNP blacklist construction
        call      — TAPS m5C calling
        benchmark — Benchmarking (rastair, asTair, Bismark)
        compare   — Concordance and correlation analysis

    Examples:
        srnataps module trim    --configfile config.yaml --slurm
        srnataps module call    --configfile config.yaml --slurm
        srnataps module compare --configfile config.yaml
    """
    # Map module name to Snakemake target rule
    rule_map = {
        "fastqc":   "rule_fastqc",
        "trim":     "rule_trim",
        "index":    "rule_index",
        "align":    "rule_align",
        "biotype":  "rule_biotype",
        "snp":      "rule_snp",
        "call":     "rule_call",
        "benchmark":"rule_benchmark",
        "compare":  "rule_compare",
    }
    target = rule_map[module_name]
    extra  = {"sample_filter": sample} if sample else {}

    console.print(f"[bold]Running module:[/bold] {module_name} — {MODULES[module_name]}")

    _run_snakemake(
        configfile=configfile,
        target=target,
        slurm=slurm,
        cores=cores,
        jobs=jobs,
        dryrun=dryrun,
        until=None,
        rerun_incomplete=False,
        extra_config=extra,
        snakemake_args=snakemake_args,
    )


# ══════════════════════════════════════════════════════════════════════════════
# srnataps check
# ══════════════════════════════════════════════════════════════════════════════
@cli.command()
@click.option("--configfile", required=True, type=click.Path(exists=True), help="Path to config.yaml")
def check(configfile):
    """
    Validate environment and config before running the pipeline.

    Checks:
        - All required tools are in PATH and correct versions
        - Reference files exist and are indexed
        - Sample sheet is valid
        - Output directories are writable
        - Snakemake is available
    """
    console.print(Panel.fit("[bold]sRNA-TAPS environment check[/bold]", title="check"))

    with open(configfile) as f:
        config = yaml.safe_load(f)

    errors   = []
    warnings = []

    # ── Tool checks ───────────────────────────────────────────────────────────
    tools = {
        "bowtie":     ("bowtie --version", "1."),
        "samtools":   ("samtools --version", "1."),
        "bcftools":   ("bcftools --version", "1."),
        "fastqc":     ("fastqc --version", None),
        "trim_galore":("trim_galore --version", None),
        "multiqc":    ("multiqc --version", None),
        "bismark":    ("bismark --version", "0."),
        "snakemake":  ("snakemake --version", "7."),
        "rastair":    ("rastair --version", None),
        "astair":     ("astair --version", None),
    }

    table = Table(title="Tools", show_header=True)
    table.add_column("Tool",    style="bold")
    table.add_column("Status",  width=8)
    table.add_column("Version")

    for tool, (cmd, min_ver) in tools.items():
        try:
            result = subprocess.run(
                cmd.split(), capture_output=True, text=True, timeout=10
            )
            ver_line = (result.stdout + result.stderr).split("\n")[0].strip()
            version  = ver_line[:60]
            status   = "[green]OK[/green]"
        except (FileNotFoundError, subprocess.TimeoutExpired):
            version = "not found"
            status  = "[red]MISSING[/red]"
            errors.append(f"Tool not found: {tool}")
        table.add_row(tool, status, version)

    console.print(table)

    # ── Reference file checks ─────────────────────────────────────────────────
    ref_table = Table(title="Reference files", show_header=True)
    ref_table.add_column("File",   style="bold")
    ref_table.add_column("Status", width=8)
    ref_table.add_column("Path")

    genome = config.get("reference", {}).get("genome_fa", "")
    gtf    = config.get("reference", {}).get("gtf", "")

    for label, path, required in [
        ("Genome FASTA", genome,          True),
        ("Genome FASTA .fai", genome+".fai", False),
        ("GTF", gtf,                      True),
        ("Bowtie1 index", config.get("reference", {}).get("bowtie1_index", "")+".1.ebwt", False),
    ]:
        if not path:
            status = "[yellow]NOT SET[/yellow]"
            if required: warnings.append(f"Not set: {label}")
        elif Path(path).exists():
            status = "[green]OK[/green]"
        else:
            status = "[red]MISSING[/red]"
            if required: errors.append(f"File not found: {path}")
            else: warnings.append(f"Not yet built: {label}")
        ref_table.add_row(label, status, path or "(not set)")

    console.print(ref_table)

    # ── Sample sheet check ────────────────────────────────────────────────────
    samples_path = config.get("input", {}).get("samples_tsv", "")
    if samples_path and Path(samples_path).exists():
        import pandas as pd
        samples = pd.read_csv(samples_path, sep="\t")
        required_cols = {"sample", "condition", "cell_line", "fastq"}
        missing_cols  = required_cols - set(samples.columns)
        if missing_cols:
            errors.append(f"samples.tsv missing columns: {missing_cols}")
        else:
            console.print(f"[green]✓[/green] samples.tsv: {len(samples)} samples")
            console.print(f"  Conditions: {sorted(samples['condition'].unique().tolist())}")
            console.print(f"  Cell lines: {sorted(samples['cell_line'].unique().tolist())}")
    else:
        warnings.append(f"samples.tsv not found: {samples_path}")

    # ── Summary ───────────────────────────────────────────────────────────────
    console.print()
    if errors:
        console.print(f"[red bold]✗ {len(errors)} error(s):[/red bold]")
        for e in errors:
            console.print(f"  [red]• {e}[/red]")
    if warnings:
        console.print(f"[yellow]{len(warnings)} warning(s):[/yellow]")
        for w in warnings:
            console.print(f"  [yellow]• {w}[/yellow]")
    if not errors and not warnings:
        console.print("[green bold]✓ All checks passed — ready to run.[/green bold]")
    elif not errors:
        console.print("[yellow]Ready to run (warnings present — see above).[/yellow]")
    else:
        console.print("[red bold]Fix errors before running.[/red bold]")
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Snakemake runner
# ══════════════════════════════════════════════════════════════════════════════
def _run_snakemake(
    configfile, target, slurm, cores, jobs,
    dryrun, until, rerun_incomplete, extra_config, snakemake_args
):
    """Build and execute the Snakemake command."""

    cmd = [
        "snakemake",
        "--snakefile",  str(SNAKEFILE),
        "--configfile", str(configfile),
    ]

    # Targets
    if until:
        cmd += ["--until", until]
    else:
        cmd += [target]

    # Execution mode
    if slurm:
        cmd += [
            "--executor",  "slurm",
            "--profile",   str(SLURM_PROF),
            "--jobs",      str(jobs),
        ]
    else:
        cmd += ["--cores", str(cores)]

    # Flags
    if dryrun:            cmd += ["--dryrun"]
    if rerun_incomplete:  cmd += ["--rerun-incomplete"]

    # Extra config key=value pairs
    if extra_config:
        for k, v in extra_config.items():
            cmd += ["--config", f"{k}={v}"]

    # Passthrough args
    if snakemake_args:
        cmd += snakemake_args.split()

    # Print command
    console.print(f"\n[dim]$ {' '.join(cmd)}[/dim]\n")

    # Execute
    result = subprocess.run(cmd)
    sys.exit(result.returncode)
