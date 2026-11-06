# =============================================================================
# rules/report.smk — R figure report generation, one rule per pipeline stage
#
# Each stage's R script reads that stage's outputs and writes figures into
# report/figures/. Because each script emits several figures (pdf/png/svg),
# each rule tracks a sentinel .done file rather than enumerating every figure.
#
# Stage → script mapping:
#   report_qc          01_qc.R           (FastQC + bowtie logs)        → rule all
#   report_biotype     02_biotype.R      (biotype summaries)          → rule all
#   report_modification 03_modification.R (TAPS call TSVs)            → rule all
#   report_seqlogos    05_sequence_logos.R (TAPS TSVs + BSgenome)     → rule all
#   report_benchmark   04_benchmark.R    (09.compare/*.tsv)           → rule all_benchmark
#   report_mirna_atlas 06_mirna_atlas.R  (stringent miRNA calls + GTF) → rule all
#   report_global_distribution 07_global_distribution.R (stringent calls) → rule all
#
# SRNATAPS_R_DIR is exported so each script can source 00_setup.R regardless
# of the working directory.
# =============================================================================

R_DIR     = Path(workflow.basedir).parent / "report" / "R"
FIG_DIR   = OUTDIR / "report" / "figures"
REPORT_LOG = LOG_DIR / "report"


rule report_qc:
    """QC figures: read-length distribution + mapping rates."""
    input:
        fastqc = expand(str(QC_DIR / "pre_trim" / "{sample}_fastqc.html"), sample=SAMPLES),
        bams   = expand(str(ALIGN_DIR / "{sample}.sorted.bam"), sample=SAMPLES),
    output:
        done = touch(str(FIG_DIR / ".report_qc.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "01_qc.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_qc.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            > {log} 2>&1
        """


rule report_biotype:
    """Biotype composition figures."""
    input:
        bams = expand(str(BIOTYPE_DIR / "{biotype}" / "{sample}_{biotype}.sorted.bam"),
                      sample=SAMPLES, biotype=BIOTYPES),
    output:
        done = touch(str(FIG_DIR / ".report_biotype.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "02_biotype.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_biotype.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            > {log} 2>&1
        """


rule report_modification:
    """TAPS modification-rate figures (distribution, top sites, context, ...)."""
    input:
        calls = expand(str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
                       sample=SAMPLES, biotype=BIOTYPES),
    output:
        done = touch(str(FIG_DIR / ".report_modification.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "03_modification.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_modification.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            > {log} 2>&1
        """


rule report_seqlogos:
    """Sequence logos around high-confidence m5C sites from the configured FASTA."""
    input:
        calls = expand(str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
                       sample=SAMPLES, biotype=BIOTYPES),
        fasta = config["reference"]["genome_fa"],
        fai = config["reference"]["genome_fa"] + ".fai",
    output:
        done = touch(str(FIG_DIR / ".report_seqlogos.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "05_sequence_logos.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_seqlogos.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            --calls-dir {CALLS_DIR} \
            --genome {input.fasta} \
            > {log} 2>&1
        """


rule report_benchmark:
    """Benchmark comparison figures (opt-in; part of all_benchmark)."""
    input:
        concordance = str(COMPARE_DIR / "concordance_summary.tsv"),
        correlation = str(COMPARE_DIR / "correlation_summary.tsv"),
    output:
        done = touch(str(FIG_DIR / ".report_benchmark.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "04_benchmark.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_benchmark.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            > {log} 2>&1
        """


rule report_mirna_atlas:
    """Complete atlas of replicate-tested miRNA modification sites, coloured by significance."""
    input:
        calls = expand(
            str(REPLICATE_DIR / "miRNA" / "treat_{cell_line}_replicate_miRNA_taps.tsv"),
            cell_line=CELL_LINES,
        ),
        gtf = config["reference"]["gtf"],
    output:
        done = touch(str(FIG_DIR / ".report_mirna_atlas.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "06_mirna_atlas.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_mirna_atlas.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            --gtf {input.gtf} \
            > {log} 2>&1
        """


rule report_global_distribution:
    """Global stage-07 modification-rate distributions by RNA biotype."""
    input:
        calls = expand(
            str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
            biotype=BIOTYPES, sample=SAMPLES,
        ),
    output:
        done = touch(str(FIG_DIR / ".report_global_distribution.done")),
    params:
        r_dir  = str(R_DIR),
        script = str(R_DIR / "07_global_distribution.R"),
        outdir = str(OUTDIR),
        figdir = str(FIG_DIR),
    log:
        str(REPORT_LOG / "report_global_distribution.log"),
    shell:
        """
        mkdir -p {params.figdir} $(dirname {log})
        SRNATAPS_R_DIR={params.r_dir} \
        Rscript {params.script} \
            --outdir {params.outdir} \
            --figdir {params.figdir} \
            --calls-dir {CALLS_DIR} \
            > {log} 2>&1
        """
