<p align="center">
  <img src="sRNA-taps_logo.png" alt="SRNA-TAPS logo">
</p>

# sRNA-TAPS
**TAPS-based m5C and 5hmC detection pipeline for small RNA sequencing**

<p align="center">
  <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  <img src="https://img.shields.io/badge/Python-3.10+-blue.svg" alt="Python">
  <img src="https://img.shields.io/badge/Platform-Linux-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/Genome-GRCh38-green.svg" alt="Genome">
  <img src="https://img.shields.io/badge/Status-In%20Development-orange.svg" alt="Status">
</p>

sRNA-TAPS pipeline is developed to detect 5-methylcytosine (m5C) and 5-hydroxymethylcytosine (5hmC) in small RNA using TET-assisted pyridine borane sequencing (TAPS). The pipeline supports miRNA, tRNA, rRNA, snoRNA, snRNA, piRNA, and lncRNA biotypes from human samples (hg38). The pipeline was adapted from the original DNA TAPS method (Liu et al., Nature Biotechnology 2019) for application to biological RNA samples without synthetic spike-in controls. TAPS operates through a two-step chemical conversion. First, ten-eleven translocation (TET) enzymes oxidise 5mC and 5hmC to 5-carboxylcytosine (5caC). Second, pyridine borane reduces 5caC to dihydrouracil (DHU), which is subsequently read as thymine during PCR amplification. The result is a C-to-T transition at modified cytosine positions in the sequencing reads. Crucially, unmodified cytosines remain unchanged and read as C. This is the inverse of bisulfite sequencing, where unmodified C is converted and modified 5mC is protected.

## Table of Contents

- [Chemistry](#chemistry)
- [Experimental Design](#experimental-design)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
- [Pipeline Steps in Detail](#pipeline-steps-in-detail)
- [The Delta (δ) Score](#the-delta-δ-score)
- [SNP Filtering](#snp-filtering)
- [Output](#output)
- [Requirements](#requirements)
- [Citation](#citation)
- [License](#license)
- [Affiliation](#affiliation)

## Chemistry

TAPS operates through a two-step chemical conversion:
1. TET enzymes oxidise m5C and 5hmC to 5caC
2. Pyridine borane reduces 5caC to dihydrouracil (DHU), read as **T** during PCR

**Unmodified C stays as C. Modified C → T.** This is the inverse of bisulfite sequencing.

### Why TAPS for small RNA?

Standard bisulfite sequencing is poorly suited for small RNA because the harsh bisulfite treatment degrades short RNA molecules, the high conversion rate (~99%) of unmodified C makes short reads difficult to align uniquely, and existing bisulfite tools assume DNA chemistry where the signal direction is inverted relative to TAPS. sRNA-TAPS uses TAPS-aware alignment parameters and a custom methylation caller that correctly interprets C-to-T transitions as genuine m5C signal.

### Why not use lambda phage spike-ins?

Lambda phage spike-ins (used in the original DNA TAPS paper) are incompatible with RNA TAPS. TET enzymes are DNA dioxygenases — their activity on RNA differs fundamentally from double-stranded DNA. RNA library preparation (adapter ligation, reverse transcription) does not process DNA spike-ins equivalently, and lambda conversion efficiency calibrates DNA chemistry, not RNA chemistry. sRNA-TAPS instead uses the `pb_Ctrl` condition for position-specific background subtraction, with known mitochondrial rRNA m5C sites serving as internal positive controls to confirm successful chemistry conversion.

## Experimental Design

sRNA-TAPS implements a three-condition control framework:

| Condition | Treatment | Purpose |
|-----------|-----------|---------|
| `treat` | TET oxidation + pyridine borane | Full TAPS — detects 5mC and 5hmC |
| `pb_Ctrl` | Pyridine borane only (no TET) | Background control — chemistry noise without TET |
| `no-treat` | No chemistry | Baseline — sequencing error rate only |

The `pb_Ctrl` samples quantify the background C-to-T conversion rate from pyridine borane chemistry alone. Subtracting this from the treated sample isolates only the TET-dependent signal, which represents genuine 5mC or 5hmC.

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

## Quick Start

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

Run individual pipeline steps independently:

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

## Pipeline Overview

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

## Pipeline Steps in Detail

### Step 1 — Quality Control

**Tool:** FastQC, MultiQC

Before any analysis, the quality of raw sequencing reads must be assessed. FastQC evaluates per-base quality scores, GC content, sequence duplication levels, overrepresented sequences, and adapter content. MultiQC aggregates results across all samples into a single interactive report. For small RNA libraries, a dominant read length peak around 18-22 nt (miRNA) and 26-32 nt (piRNA/tRNA fragments) is expected, along with adapter sequences since small RNA inserts are shorter than the sequencing read length.

### Step 2 — Adapter Trimming

**Tool:** TrimGalore (wrapper for Cutadapt)

Small RNA sequencing libraries are prepared by ligating adapters to the 3' end of RNA molecules. Because miRNAs are 18-22 nt and the sequencing read length is typically 50-75 nt, every read contains adapter sequence after the insert. If adapters are not removed they prevent alignment or cause misalignment to the genome. TrimGalore automatically detects adapter sequences and applies quality trimming simultaneously. A minimum length filter of 15 nt is applied to remove reads too short to align reliably.

### Step 3 — TAPS-aware Alignment

**Tool:** Bowtie 1.3.1, SAMtools

TAPS data requires careful alignment parameter choices. The alignment must be configured to tolerate C-to-T mismatches in the reads, because these transitions represent genuine 5mC modifications rather than sequencing errors or alignment artefacts. Four critical parameters are used:

| Parameter | Purpose |
|-----------|---------|
| `-v 2` | Allows up to 2 total mismatches without quality weighting, so methylated C→T transitions are tolerated rather than penalised |
| `--norc` | Restricts alignment to the forward strand only — small RNA libraries are strand-specific, and reverse complement alignment would introduce false G→A calls incorrectly interpreted as modified cytosines |
| `-k 10 --best --strata` | Reports up to 10 alignments per read from the best-scoring stratum, writing the `XA:i:N` tag recording multi-mapping count — essential for fractional weighting in the methylation caller |
| `--sam` | Required for the XA tag to be written |

Alignment is performed to GRCh38 at the genome level rather than a transcriptome, because small RNA species including piRNAs, tRNAs (~600 gene copies), and snoRNAs are encoded at specific genomic loci whose coordinates are needed for downstream annotation.

### Step 4 — Biotype Annotation and BAM Splitting

**Tool:** Custom Python script using pysam, Ensembl GRCh38 v112 GTF

Different RNA classes are biologically distinct and require separate interpretation. An m5C site in a tRNA at a known structural position means something completely different from an m5C site in a miRNA seed region. Mixing all RNA types in a single methylation call file makes biological interpretation impossible. Each aligned read is intersected with Ensembl gene annotations and assigned to the highest-priority biotype following a strict hierarchy:

```
miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other
```

This priority order ensures reads overlapping multiple annotations are assigned to the most biologically specific category.

> **Note:** TAPS chemistry is known to alter small RNA library composition relative to untreated samples, particularly enriching for rRNA. Users should expect biotype proportions to differ between treated and untreated conditions. This is a technical consequence of the pyridine borane step and should be accounted for in downstream interpretation.

### Step 5 — TAPS Methylation Calling

**Tool:** Custom Python caller (`taps_calling_fast.py`) using pysam, multiprocessing

The methylation caller implements TAPS logic at single-base resolution. For each cytosine position covered by at least one read, it counts reads showing C (unmodified) and T (modified after TAPS conversion):

```
mod_rate = T_count / (C_count + T_count)
```

Key design decisions that make this caller correct for TAPS RNA data:

| Feature | Implementation |
|---------|---------------|
| **Strand-specific logic** | Forward strand: modified C reads as T. Reverse strand: modified C appears as G in reference and reads as A after conversion. Both orientations handled separately. |
| **CIGAR-aware parsing** | Uses `get_aligned_pairs()` from pysam to correctly handle insertions, deletions, and soft-clipped bases, preventing position misalignment |
| **Multi-mapper fractional weighting** | Reads with `XA:i:N` tag contribute `1/N` to each location rather than 1.0, preventing inflation at multi-mapping positions such as tRNA gene copies and rRNA repeats |
| **Base quality filtering** | Bases with Phred quality < 20 are excluded, preventing low-quality calls from contributing false C→T transitions |
| **Chromosomal parallelisation** | Work is distributed across chromosomes using Python multiprocessing, reducing wall time proportionally to available CPU cores |

Minimum coverage thresholds per biotype: **rRNA 10x**, **miRNA and snoRNA 5x**, **tRNA 3x**.

Output columns: `chrom, start, end, context, mod_count, unmod_count, coverage, mod_rate`

### Step 6 — Background Correction and Replicate Merging

**Tool:** Custom Python script

The `pb_Ctrl` samples define the chemistry background rate at each genomic position. Background correction is performed site by site:

```
δ = treat_mod_rate − pb_Ctrl_mod_rate
```

This subtraction removes two sources of false signal: the intrinsic rate of pyridine borane-mediated C→T conversion at unmodified cytosines, and position-specific biases in the chemistry related to local sequence context. Sites detected in fewer than two replicates are discarded — a modification detectable in only one replicate cannot be distinguished from sample-specific noise. This replicate reproducibility filter is the most important quality gate in the pipeline.

### Step 7 — Differential Methylation Analysis

**Tool:** Custom Python cross-condition comparison

Genuine m5C sites are identified using three simultaneous criteria:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| `delta > 0.1` | δ > 10% | Accounts for residual chemistry noise; ensures only meaningful signal is reported |
| `notx_mean < 0.05` | no-treat < 5% | Excludes C→T SNPs and RNA editing sites, which show signal in all conditions regardless of TAPS chemistry |
| `rep >= 2` | ≥ 2/3 replicates | Strongest filter against false positives; stochastic noise is not reproducible across replicates |

Sites passing all three criteria are classified by confidence:

| Confidence | Criteria |
|------------|----------|
| **High** | rep = 3/3 and δ > 0.3 |
| **Medium** | rep ≥ 2 and δ > 0.15 |
| **Low** | rep = 2 and δ > 0.1 |

### Step 8 — Genomic Annotation

**Tool:** bedtools intersect, Ensembl GRCh38 v112 GTF

Genomic coordinates alone are biologically uninterpretable. Each candidate m5C site is intersected with the Ensembl gene annotation to assign gene name, gene biotype, and feature type. For miRNA candidates, gene names follow miRBase nomenclature. For sites overlapping multiple gene annotations, the most specific annotation is retained using the priority hierarchy from Step 4.

## The Delta (δ) Score

### What is δ?

Delta is the background-corrected methylation signal — how much of the C-to-T conversion observed in the treated sample is genuinely due to 5mC, after subtracting the noise introduced by the chemistry itself.

```
δ = treat_mod_rate − pb_Ctrl_mod_rate
```

### Why subtract pb_Ctrl?

Pyridine borane does not act exclusively on 5-carboxylcytosine. It causes a low level of non-specific C-to-T conversion at unmodified cytosines (~8-9% background rate). The pb_Ctrl samples went through pyridine borane but not TET oxidation, meaning any C-to-T signal in pb_Ctrl is pure chemistry noise. Subtracting it isolates only the TET-dependent signal representing genuine 5mC or 5hmC.

### Concrete example

```
treat_mod_rate   = 0.823  (82.3% of reads show T at this position)
pb_Ctrl_mod_rate = 0.323  (32.3% background from pyridine borane alone)

δ = 0.823 − 0.323 = 0.500
```

Without subtraction, one might conclude 82% methylation. But 32% is chemistry noise. The genuine 5mC fraction is 50%.

### Interpreting δ values

| δ value | Interpretation |
|---------|---------------|
| ~0.00 | No methylation above background — unmodified cytosine |
| 0.10–0.20 | Low-level methylation — weak but reproducible signal |
| 0.20–0.40 | Moderate methylation — confident candidate |
| 0.40–0.70 | High methylation — strong m5C signal |
| > 0.70 | Very high or near-stoichiometric methylation |

### Important limitation

δ is a relative measure, not absolute stoichiometry. Without a fully methylated RNA spike-in control, conversion efficiency cannot be determined. δ values are internally consistent and sufficient for site identification and cell-line comparison, but should not be interpreted as absolute methylation fractions.

## SNP Filtering

sRNA-TAPS applies three-layer polymorphism filtering **before** counting C→T events. This is critical: a C/T heterozygous SNP is chemically indistinguishable from a TAPS m5C signal.

| Layer | Method | Flag |
|-------|--------|------|
| 1 | dbSNP common C→T/G→A variants (AF ≥ 1%) | `SNP_KNOWN` |
| 2 | Cell-line-specific SNPs from no-treat BAMs | `SNP_SAMPLE` |
| 3 | No-treat C→T rate ≥ 40% (heterozygosity) | `SNP_HET` |

The `snp_flag` column in output TSVs records the flag for every site. Only `PASS` sites proceed to statistical testing.

Additionally, sRNA-TAPS supports cross-validation with **Rastair v2.1.1**, which applies an independent machine-learning-based SNP correction at each called position.

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

MIT License — see [LICENSE](LICENSE) for details.

Copyright (c) 2026 Bennett Henzeler, Institute of Chemical Epigenetics, Ludwig-Maximilians-Universität München

## Affiliation

<p align="left">
  <img src="Lab_Logo.png" height="80" alt="Institute of Chemical Epigenetics-Munich">
  <br><br>
  <a href="https://schneider.cup.uni-muenchen.de">Schneider Lab</a><br>
  Institute of Chemical Epigenetics-Munich (ICEM)<br>
  Ludwig-Maximilians-Universität München<br>
  Munich, Germany
</p>
