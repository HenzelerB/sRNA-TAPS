# =============================================================================
# rules/report.smk — MultiQC and final HTML report
# =============================================================================

rule multiqc:
    """Aggregate all QC metrics with MultiQC."""
    input:
        expand(str(QC_DIR / "pre_trim"  / "{sample}_fastqc.html"), sample=SAMPLE_NAMES),
        expand(str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.html"), sample=SAMPLE_NAMES),
        expand(str(ALIGN_DIR / "{sample}.sorted.bam"), sample=SAMPLE_NAMES),
    output:
        html = str(REPORT_DIR / "multiqc_report.html"),
    params:
        search_dirs = f"{QC_DIR} {TRIM_DIR} {ALIGN_DIR}",
        outdir      = str(REPORT_DIR),
        title       = config["project"]["name"],
    log:
        str(LOG_DIR / "report" / "multiqc.log"),
    resources:
        mem_mb  = 8000,
        runtime = 60,
    shell:
        """
        mkdir -p {params.outdir}
        multiqc \
            {params.search_dirs} \
            --outdir {params.outdir} \
            --filename multiqc_report \
            --title "{params.title} — sRNA-TAPS QC" \
            --comment "TAPS small RNA methylation pipeline" \
            --force \
            > {log} 2>&1
        """
