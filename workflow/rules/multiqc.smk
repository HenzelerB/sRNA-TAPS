# =============================================================================
# rules/multiqc.smk — MultiQC aggregation report
# =============================================================================

rule multiqc:
    """Aggregate FastQC, Trim Galore, and alignment QC into a single MultiQC report."""
    input:
        expand(str(QC_DIR / "pre_trim"  / "{sample}_fastqc.html"), sample=SAMPLES),
        expand(str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.html"), sample=SAMPLES),
        expand(str(ALIGN_DIR / "{sample}.sorted.bam"), sample=SAMPLES),
    output:
        html = str(OUTDIR / "report" / "multiqc_report.html"),
    params:
        outdir   = str(OUTDIR / "report"),
        trim_dir = str(TRIM_DIR),
        qc_dir   = str(QC_DIR),
        log_dir  = str(LOG_DIR),
    log:
        str(LOG_DIR / "multiqc.log"),
    shell:
        """
        mkdir -p {params.outdir}
        multiqc \
            {params.qc_dir} \
            {params.trim_dir} \
            {params.log_dir} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            --title "sRNA-TAPS QC report" \
            > {log} 2>&1
        """
