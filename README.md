<p align="center">
  <img src="sRNA-taps_logo.png" alt="sRNA-TAPS logo" width="700">
</p>

<h1 align="center">sRNA-TAPS</h1>
<p align="center"><strong>TAPS-based m5C and 5hmC detection pipeline for small RNA sequencing</strong></p>

<p align="center">
  <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  <img src="https://img.shields.io/badge/Python-3.10+-blue.svg" alt="Python">
  <img src="https://img.shields.io/badge/Platform-Linux-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/Genome-GRCh38-green.svg" alt="Genome">
  <img src="https://img.shields.io/badge/Status-In%20Development-orange.svg" alt="Status">
  <img src="https://img.shields.io/badge/Tests-66%2F66%20passing-brightgreen.svg" alt="Tests">
</p>

<p align="center">
sRNA-TAPS detects 5-methylcytosine (m5C) and 5-hydroxymethylcytosine (5hmC) in small RNA using TET-assisted pyridine borane sequencing (TAPS). It covers miRNA, tRNA, rRNA, snoRNA, snRNA, piRNA, and lncRNA biotypes from human samples (hg38), and was developed from the original DNA TAPS method (Liu et al., <em>Nature Biotechnology</em> 2019) for biological RNA without synthetic spike-in controls.
</p>

---

## 📋 Table of Contents

<ol>
  <li><a href="#-chemistry">Chemistry</a></li>
  <li><a href="#-experimental-design">Experimental Design</a></li>
  <li><a href="#️-installation">Installation</a></li>
  <li><a href="#-quick-start">Quick Start</a></li>
  <li><a href="#-test-dataset">Test Dataset</a></li>
  <li><a href="#-usage">Usage</a></li>
  <li><a href="#-pipeline-overview">Pipeline Overview</a></li>
  <li><a href="#-pipeline-steps-in-detail">Pipeline Steps in Detail</a></li>
  <li><a href="#-the-delta-δ-score">The Delta (δ) Score</a></li>
  <li><a href="#️-snp-filtering">SNP Filtering</a></li>
  <li><a href="#-output">Output</a></li>
  <li><a href="#-requirements">Requirements</a></li>
  <li><a href="#-citation">Citation</a></li>
  <li><a href="#️-license">License</a></li>
  <li><a href="#️-affiliation">Affiliation</a></li>
</ol>

---

## 🧪 Chemistry

TAPS works through two sequential reactions: TET enzymes oxidise m5C and 5hmC to 5-carboxylcytosine (5caC), which pyridine borane then reduces to dihydrouracil (DHU). DHU is read as T during PCR. Unmodified cytosines pass through the chemistry unchanged, so the readout is the opposite of bisulfite sequencing — a C→T transition marks a modified base rather than an unmodified one.

#### Why TAPS for small RNA?

Bisulfite treatment degrades short RNA molecules and converts ~99% of unmodified cytosines, making 18–22 nt reads nearly impossible to align uniquely. TAPS leaves unmodified cytosines intact, preserving alignment specificity. The pipeline uses Bowtie1 parameters tuned to tolerate the C→T mismatches introduced by modification rather than penalising them as sequencing errors.

#### Why not use lambda phage spike-ins?

Lambda spike-ins calibrate conversion efficiency for double-stranded DNA. TET enzymes act differently on RNA, and the adapter ligation and reverse transcription steps in small RNA library preparation do not process DNA spike-ins the same way. sRNA-TAPS uses `pb_Ctrl` samples for position-specific background correction instead, with mitochondrial rRNA m5C sites (C1402, C1407 in 12S rRNA) as internal positive controls for chemistry conversion.

---

## 🔬 Experimental Design

Three conditions are required:

| Condition | Treatment | Purpose |
|-----------|-----------|---------|
| `treat` | TET oxidation + pyridine borane | Full TAPS — detects 5mC and 5hmC |
| `pb_Ctrl` | Pyridine borane only (no TET) | Captures chemistry background without TET |
| `no-treat` | No chemistry | Sequencing error baseline |

The `pb_Ctrl` condition is essential. Pyridine borane converts a small fraction of unmodified cytosines non-specifically (~8–9%), and this background varies by sequence context. Subtracting `pb_Ctrl` rates from `treat` at each position leaves only the TET-dependent signal.

---

## ⚙️ Installation

#### Conda (recommended)

```bash
conda install -c bioconda -c conda-forge srna-taps
```

#### pip

```bash
pip install sRNA-TAPS
```

#### From source

```bash
git clone https://github.com/HenzelerB/sRNA-TAPS
cd sRNA-TAPS
pip install -e .
```

#### Full environment

```bash
conda env create -f environment.yaml
conda activate sRNA-TAPS
```

#### R packages (for report figures)

```bash
Rscript install_R_packages.R
```

---

## 🚀 Quick Start

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

---

## 🧬 Test Dataset

The repository includes a synthetic FASTQ simulator (`tests/simulate_taps_srna.py`) that generates small RNA reads with realistic TAPS chemistry. Known m5C positions are seeded into the templates at 40–80% conversion in the `treat` condition and 2–5% background in `no_treat`, so you can verify end-to-end pipeline behaviour without real sequencing data.

#### Generating test FASTQs

```bash
python3 tests/simulate_taps_srna.py \
    --outdir /path/to/test_fastq \
    --reads  100000 \
    --seed   42
```

This writes **9 samples** (3 conditions × 3 replicates, HEK cell line):

| Sample | Condition | Description |
|--------|-----------|-------------|
| `no-treat_Ctrl_HEK_R1/R2/R3` | `no_treat` | No chemistry — sequencing error baseline |
| `pb_Ctrl_HEK_R1/R2/R3` | `pb_ctrl` | PB only — chemistry background without TET |
| `treat_HEK_R1/R2/R3` | `treat` | TET + PB — genuine TAPS signal |

Reads are 18–50 nt with a TruSeq small RNA 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`). Biotype proportions: miRNA 40%, tRNA 25%, rRNA 25%, snoRNA 10%, drawn from real hg38 loci (hsa-miR-21-5p, mt-tRNA-Leu, mt-12S rRNA, SNORD14). A `samples.tsv` is written alongside the FASTQs.

#### Running the test pipeline

These files use standard hg38 coordinates and work directly with an existing Bowtie1 index and Ensembl GTF. Start from trimming:

```bash
# Step 02: Adapter trimming
trim_galore --small_rna --length 18 --max_length 50 \
    --output_dir 03.trimGalore /path/to/test_fastq/*.fastq.gz

# Step 03: Bowtie1 alignment
bowtie -x /path/to/hg38_index -q sample_trimmed.fq.gz \
    --norc -v 2 -k 10 --best --strata -m 100 -S sample.sam

# Step 04: Biotype annotation
python3 05_annotate_biotype.py \
    --bam sample.sorted.bam --gtf hg38.gtf \
    --out_dir 05.biotype_bams --sample <sample_name>

# Step 05: Cell-line SNP blacklist
samtools merge notreated_HEK_merged.bam \
    no-treat_Ctrl_HEK_R{1,2,3}.sorted.bam
python3 06_build_snp_blacklist.py \
    --bam notreated_HEK_merged.bam --fasta hg38.fa \
    --out 06.snp_resources/sample_snps_HEK.bed \
    --min-af 0.20 --min-cov 5 --cell-line HEK

# Step 06: TAPS modification calling
python3 07_taps_calling.py \
    --bam 05.biotype_bams/miRNA/treat_HEK_R1_miRNA.sorted.bam \
    --fasta hg38.fa \
    --out 07.taps_calls/miRNA/treat_HEK_R1_miRNA_taps.tsv \
    --min-cov 3 --context ALL \
    --sample-snp-bed 06.snp_resources/sample_snps_HEK.bed \
    --cell-line HEK
```

#### Expected results

Around 80% of trimmed reads align to hg38. At the two seeded miRNA m5C sites:

| Position | Gene | Condition | mod_rate |
|----------|------|-----------|----------|
| chr17:59841313 | hsa-miR-21-5p | TET+PB | **90.9%** |
| chr17:59841313 | hsa-miR-21-5p | Untreated | **1.7%** |
| chrX:66018955 | hsa-miR-223-3p | TET+PB | **71.6%** |
| chrX:66018955 | hsa-miR-223-3p | Untreated | **0.2%** |

```bash
# Quick check — expect ~0.9 in treat, ~0.02 in no_treat
grep "^17.*59841313" 07.taps_calls/miRNA/treat_HEK_R1_miRNA_taps.tsv
grep "^17.*59841313" 07.taps_calls/miRNA/no-treat_Ctrl_HEK_R1_miRNA_taps.tsv
```

> **Multi-lane data:** If your sequencing run was split across lanes, merge the per-lane FASTQs before running: `cat sample_L001.fastq.gz sample_L002.fastq.gz > sample_merged.fastq.gz`. The test dataset is already single-file per sample.

---

## 📖 Usage

```
Usage: srnataps [OPTIONS] COMMAND [ARGS]...

  sRNA-TAPS: TAPS-based m5C detection for small RNA sequencing.

Commands:
  init    Initialise a new project (config.yaml + samples.tsv)
  run     Run the full pipeline
  module  Run a single pipeline module
  check   Validate environment and config
```

#### `srnataps run`

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

#### `srnataps module`

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

---

## 🔄 Pipeline Overview

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
09. Report          R figures (PDF/PNG/SVG) · interactive HTML (Plotly)
```
*requires `--benchmark` flag

---

## 🔍 Pipeline Steps in Detail

### Step 1 — Quality Control
**Tool:** FastQC, MultiQC

FastQC is run on raw FASTQs and MultiQC aggregates the results across all samples. For small RNA libraries, expect a dominant length peak at 18–22 nt (miRNA) and a second peak at 26–32 nt (tRNA fragments and piRNAs). All reads will contain adapter sequence, since small RNA inserts are shorter than the read length.

---

### Step 2 — Adapter Trimming
**Tool:** TrimGalore (Cutadapt)

TrimGalore removes the TruSeq small RNA 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`) and quality-trims the 3′ end. Reads shorter than 18 nt or longer than 50 nt are discarded. The `--small_rna` flag sets sensible defaults for this library type.

---

### Step 3 — TAPS-aware Alignment
**Tool:** Bowtie 1.3.1, SAMtools

Alignment is to the whole GRCh38 genome rather than a transcriptome — tRNA families (~600 gene copies), piRNA clusters, and snoRNA loci need genomic coordinates for correct downstream annotation. Four parameters are set specifically for TAPS data:

| Parameter | Purpose |
|-----------|---------|
| `-v 2` | Up to 2 total mismatches — tolerates C→T modifications without penalising them |
| `--norc` | Forward strand only — reverse complement alignment would introduce false G→A calls |
| `-k 10 --best --strata` | Up to 10 alignments per read from the best stratum; writes `XA:i:N` tag used for fractional weighting |
| `--sam` | Required for XA tag output |

---

### Step 4 — Biotype Annotation and BAM Splitting
**Tool:** Custom Python script (pysam), Ensembl GRCh38 v112 GTF

Reads are assigned to the highest-priority overlapping biotype:

```
miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other
```

Keeping biotypes in separate BAMs matters because m5C at a tRNA wobble position has completely different biological meaning from m5C in a miRNA seed region, and coverage thresholds appropriate for one biotype are not appropriate for another.

> **Note:** TAPS enriches for rRNA relative to untreated libraries due to the pyridine borane step. Biotype proportions will differ between conditions — this is expected and not a sign of a failed experiment.

---

### Step 5 — TAPS Methylation Calling
**Tool:** Custom Python caller (pysam, multiprocessing)

For each cytosine position with sufficient coverage:

```
mod_rate = T_count / (C_count + T_count)
```

| Feature | Implementation |
|---------|---------------|
| **Strand-specific logic** | Forward strand: modified C reads as T. Reverse strand: modified C reads as A. |
| **CIGAR-aware parsing** | `get_aligned_pairs()` handles indels and soft-clipping without positional offset errors |
| **Multi-mapper weighting** | Reads with `XA:i:N` contribute `1/N` per locus rather than 1, preventing inflation at tRNA gene copies and rRNA repeats |
| **Base quality filter** | Bases with Phred Q < 20 are excluded |
| **Parallel processing** | Jobs distributed across chromosomes |
| **Statistical testing** | Binomial test against background rate, BH FDR correction |

Output columns: `chrom, start, end, context, mod_count, unmod_count, coverage, mod_rate, pvalue, padj, snp_flag`

---

### Step 6 — Background Correction and Replicate Merging
**Tool:** Custom Python script

```
δ = treat_mod_rate − pb_Ctrl_mod_rate
```

Subtracting `pb_Ctrl` removes both the global pyridine borane background and sequence-context-specific biases. Sites seen in only one replicate are discarded — stochastic noise is not reproducible, so this filter is the most effective way to reduce false positives.

---

### Step 7 — Differential Methylation Analysis
**Tool:** Custom Python cross-condition comparison

A site is called as genuine m5C if it passes all three filters simultaneously:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| `delta > 0.1` | δ > 10% | Removes residual chemistry noise after background subtraction |
| `notx_mean < 0.05` | no-treat < 5% | Removes SNPs and A-to-I editing sites, which are condition-independent |
| `rep >= 2` | ≥ 2/3 replicates | Reproducibility across replicates is the strongest indicator of genuine signal |

Called sites are then stratified by confidence:

| Confidence | Criteria |
|------------|----------|
| **High** | rep = 3/3 and δ > 0.3 |
| **Medium** | rep ≥ 2 and δ > 0.15 |
| **Low** | rep = 2 and δ > 0.1 |

---

### Step 8 — Genomic Annotation
**Tool:** bedtools intersect, Ensembl GRCh38 v112 GTF

Called sites are intersected with the Ensembl annotation to assign gene name, biotype, and feature. miRNA gene names follow miRBase nomenclature. Where a site overlaps multiple features, the most specific annotation is kept using the priority hierarchy from Step 4.

---

### Step 9 — Report Generation
**Tools:** R (ggplot2, ggseqlogo, BSgenome.Hsapiens.UCSC.hg38), Python (Plotly)

The pipeline generates publication-ready static figures (PDF/PNG/SVG at 300 dpi, Arial 8pt) and a self-contained interactive HTML report. Figures cover QC, biotype composition, modification rates, condition comparisons, benchmarking concordance, and sequence logos at ±5 and ±10 nt windows around called sites.

```bash
RDIR=/path/to/sRNA-TAPS/srnataps/report/R
export SRNATAPS_R_DIR=$RDIR

# Full report
Rscript $RDIR/run_all.R \
    --outdir  /path/to/project \
    --figdir  /path/to/project/report/figures \
    --scripts $RDIR

# Skip sections you don't need
Rscript $RDIR/run_all.R --outdir /path/to/project --scripts $RDIR \
    --skip-qc --skip-bio --skip-mod --skip-bench --skip-logos

# Interactive HTML
python3 srnataps/report.py \
    --outdir /path/to/project \
    --out    /path/to/project/report/srnataps_report.html
```

---

## Δ The Delta (δ) Score

δ is the background-corrected modification rate — the fraction of C→T signal in the treated sample that is attributable to TET-dependent oxidation rather than pyridine borane chemistry alone:

```
δ = treat_mod_rate − pb_Ctrl_mod_rate
```

For example, if a position shows 82.3% C→T in `treat` and 32.3% in `pb_Ctrl`, then δ = 0.50 — half the signal is genuine 5mC and half is chemistry noise. Without this correction, the site would appear 82% methylated.

#### Interpreting δ values

| δ value | Interpretation |
|---------|---------------|
| ~0.00 | No methylation above background |
| 0.10–0.20 | Low-level methylation |
| 0.20–0.40 | Moderate methylation — confident candidate |
| 0.40–0.70 | High methylation |
| > 0.70 | Near-stoichiometric methylation |

> δ is not absolute stoichiometry. Without a fully modified RNA standard, the actual conversion efficiency is unknown. Use δ for site discovery and cross-condition comparisons, not as an absolute measure of methylation fraction.

---

## 🛡️ SNP Filtering

A C/T heterozygous SNP produces exactly the same signal as a TAPS m5C site. sRNA-TAPS filters at three levels before any C→T counts are made:

| Layer | Method | Flag |
|-------|--------|------|
| 1 | dbSNP common C→T/G→A variants (AF ≥ 1%) | `SNP_KNOWN` |
| 2 | Cell-line-specific C→T variants called from no-treat BAMs | `SNP_SAMPLE` |
| 3 | No-treat C→T rate ≥ 40% at a position (heterozygosity proxy) | `SNP_HET` |

Filtering happens upstream of the binomial test — SNP-flagged positions are excluded entirely rather than called and then filtered. The `snp_flag` column in every output TSV records the reason for exclusion. Only `PASS` sites are tested. Cross-validation with Rastair v2.1.1 adds an independent machine-learning-based SNP correction layer for benchmarking purposes.

---

## 📁 Output

```
outdir/
├── 02.fastqc/          FastQC HTML reports (pre-trim)
├── 03.trimGalore/      Trimmed FASTQs + trim reports
├── 04a.genome/         Bowtie1 genome index
├── 04b.aligned/        Sorted BAMs per sample
├── 05.biotype_bams/    Per-biotype BAMs + biotype_composition_all_samples.tsv
├── 06.snp_resources/   SNP blacklists per cell line
├── 07.taps_calls/      Per-biotype TAPS TSVs
│                       (chrom, start, end, context, mod_count, unmod_count,
│                        coverage, mod_rate, pvalue, padj, snp_flag)
├── 08.benchmark/       rastair · asTair · Bismark outputs
├── 09.compare/         concordance_summary.tsv · correlation_summary.tsv
└── report/
    ├── figures/        PDF + PNG + SVG per figure
    │   ├── 01_qc       Read length distribution, mapping rates
    │   ├── 02_bio      Biotype composition
    │   ├── 03_mod      Modification rates, top sites, condition comparison,
    │   │               waterfall, trinucleotide context
    │   ├── 04_bench    Concordance heatmap, correlation, site overlap
    │   └── 05_logos    Sequence logos (±5 and ±10 nt, per biotype)
    └── srnataps_report.html   Interactive HTML report (Plotly)
```

---

## 📦 Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥ 3.10 | Core pipeline |
| bowtie | 1.3.1 | TAPS-aware alignment |
| samtools | ≥ 1.20 | BAM processing |
| bcftools | ≥ 1.20 | Variant calling |
| fastqc | ≥ 0.12 | Quality control |
| trim-galore | ≥ 0.6.10 | Adapter trimming |
| multiqc | ≥ 1.21 | QC aggregation |
| snakemake | ≥ 7.0 | Workflow management |
| R | ≥ 4.3 | Report figures |
| bismark | ≥ 0.24 | Benchmarking only |
| rastair | ≥ 2.1 | Benchmarking only |
| asTair | ≥ 3.3 | Benchmarking only |

---

## 📄 Citation

If you use sRNA-TAPS, please cite:

> Henzeler B. et al. sRNA-TAPS: TAPS-based m5C detection in small RNA. *(in preparation)*

And the underlying TAPS method:

> Liu Y. et al. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution. *Nature Biotechnology* 37, 424–429 (2019). https://doi.org/10.1038/s41587-019-0041-2

---

## ⚖️ License

MIT License — see [LICENSE](LICENSE) for details.

Copyright (c) 2026 Bennett Henzeler, Institute of Chemical Epigenetics, Ludwig-Maximilians-Universität München

---

## 🏛️ Affiliation

<p align="left">
  <img src="Lab_Logo.png" height="80" alt="Institute of Chemical Epigenetics-Munich">
  <br><br>
  <a href="https://schneider.cup.uni-muenchen.de">Schneider Lab</a><br>
  Institute of Chemical Epigenetics-Munich (ICEM)<br>
  Ludwig-Maximilians-Universität München<br>
  Munich, Germany
</p>
