# sRNA-TAPS
<p align="center">
  <img src="sRNA-taps_logo.tif" alt="SRNA-TAPS logo">
</p>

**TAPS-based m5C and 5hmC detection pipeline for small RNA sequencing**

sRNA-TAPS detects 5-methylcytosine (m5C) and 5-hydroxymethylcytosine (5hmC) in small RNA using TET-assisted pyridine borane sequencing (TAPS). The pipeline supports miRNA, tRNA, rRNA, snoRNA, snRNA, piRNA, and lncRNA biotypes from human samples (hg38).

## Chemistry

TAPS operates through a two-step chemical conversion:
1. TET enzymes oxidise m5C and 5hmC to 5-carboxylcytosine (5caC)
2. Pyridine borane reduces 5caC to dihydrouracil (DHU), read as **T** during PCR

**Unmodified C stays as C. Modified C → T.** This is the inverse of bisulfite sequencing.

## Installation

### Conda (recommended)

```bash
conda install -c bioconda -c conda-forge srna-taps
```

### pip

```bash
pip install sRNA-TAPS
```

### From source

```bash
git clone https://github.com/bhenzeler/sRNA-TAPS
cd sRNA-TAPS
pip install -e .
```

### Full environment

```bash
conda env create -f environment.yaml
conda activate sRNA-TAPS
```

## Quick start

```bash
# 1. Initialise a new project
srnataps init \
    --outdir  ~/my_taps_project \
    --genome  /path/to/hg38.fa \
    --gtf     /path/to/hg38.gtf

# 2. Edit the sample sheet
nano ~/my_taps_project/samples.tsv

# 3. Validate environment
srnataps check --configfile ~/my_taps_project/config.yaml

# 4. Run full pipeline on SLURM cluster
srnataps run --configfile ~/my_taps_project/config.yaml --slurm

# 5. Run with benchmarking
srnataps run --configfile ~/my_taps_project/config.yaml --slurm --benchmark
```

## Usage

```
Usage: srnataps [OPTIONS] COMMAND [ARGS]...

  sRNA-TAPS: TAPS-based m5C detection for small RNA sequencing.

Commands:
  init    Initialise a new project (config.yaml + samples.tsv)
  run     Run the full pipeline
  module  Run a single pipeline module
  check   Validate environment and config
```

### `srnataps run`

```
Options:
  --configfile  PATH   Path to config.yaml  [required]
  --slurm              Submit jobs to SLURM cluster
  --benchmark          Also run rastair, asTair, and Bismark benchmarking
  --cores       INT    Local cores (ignored if --slurm)
  --jobs        INT    Max concurrent SLURM jobs [default: 50]
  --dryrun, -n         Show what would run without executing
  --until       RULE   Run up to and including this rule
```

### `srnataps module`

```
Available modules:
  fastqc    FastQC on raw merged FASTQs
  trim      Trim Galore adapter trimming (TruSeq small RNA)
  index     Bowtie1 genome index build
  align     Bowtie1 alignment
  biotype   RNA biotype BAM splitting
  snp       SNP blacklist construction
  call      TAPS m5C modification calling
  benchmark Benchmarking (rastair, asTair, Bismark)
  compare   Concordance and correlation analysis

Example:
  srnataps module call --configfile config.yaml --slurm
```

## Pipeline overview

```
rawfiles/           Raw merged FASTQs (SE, TruSeq small RNA)
    ↓
01. FastQC          Pre-trim QC
02. Trim Galore     TruSeq SR adapter trimming (--small_rna)
03. Bowtie1         Alignment: -v2 --norc -k10 --best --strata -m100
04. Biotype split   miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other
05. SNP filter      3-layer: dbSNP + cell-line-specific + heterozygosity
06. TAPS calling    Pileup → C→T counting → binomial test → BH FDR
    ↓
07. Benchmark*      rastair (Bowtie1 BAMs) · asTair (biotype BAMs) · Bismark
08. Compare*        Concordance · Pearson/Spearman correlation
09. Report          MultiQC · per-biotype TSVs · HTML summary
```
*requires `--benchmark` flag

## SNP filtering

sRNA-TAPS applies three-layer polymorphism filtering **before** counting C→T events. This is critical: a C/T heterozygous SNP is chemically indistinguishable from a TAPS m5C signal.

| Layer | Method | Flag |
|-------|--------|------|
| 1 | dbSNP common C→T/G→A variants (AF ≥ 1%) | `SNP_KNOWN` |
| 2 | Cell-line-specific SNPs from no-treat BAMs | `SNP_SAMPLE` |
| 3 | No-treat C→T rate ≥ 40% (heterozygosity) | `SNP_HET` |

The `snp_flag` column in output TSVs records the flag for every site. Only `PASS` sites proceed to statistical testing.

## Output

```
outdir/
├── 02.fastqc/          FastQC HTML reports (pre + post trim)
├── 03.trimGalore/      Trimmed FASTQs + trim reports
├── 04a.genome/         Bowtie1 genome index
├── 04b.aligned/        Sorted BAMs per sample
├── 05.biotype_bams/    Per-biotype BAMs + composition summary
├── 06.snp_resources/   SNP blacklists per cell line
├── 07.taps_calls/      Per-biotype TAPS TSVs (chrom,start,end,context,
│                       mod_count,unmod_count,coverage,mod_rate,pvalue,padj,snp_flag)
├── 08.benchmark/       rastair · asTair · Bismark outputs
├── 09.compare/         concordance_summary.tsv · correlation_summary.tsv
└── report/             multiqc_report.html
```

## Requirements

- Python ≥ 3.10
- bowtie 1.3.1
- samtools ≥ 1.20
- bcftools ≥ 1.20
- fastqc ≥ 0.12
- trim-galore ≥ 0.6.10
- multiqc ≥ 1.21
- bismark ≥ 0.24 (benchmarking)
- rastair ≥ 2.1 (benchmarking)
- asTair ≥ 3.3 (benchmarking)
- snakemake ≥ 7.0

## Citation

If you use sRNA-TAPS, please cite:

> Henzeler B. et al. sRNA-TAPS: TAPS-based m5C detection in small RNA. (in preparation)

And the underlying TAPS method:

> Liu Y. et al. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution. *Nature Biotechnology* 37, 424–429 (2019). https://doi.org/10.1038/s41587-019-0041-2

## License

MIT
