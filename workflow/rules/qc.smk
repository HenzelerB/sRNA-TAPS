# =============================================================================
# rules/qc.smk — FastQC pre- and post-trim
# =============================================================================

rule fastqc_pretrim:
    input:
        fq = lambda wc: SAMPLES_DF.loc[wc.sample, "fastq"],
    output:
        html = str(QC_DIR / "pre_trim" / "{sample}_fastqc.html"),
        zip  = str(QC_DIR / "pre_trim" / "{sample}_fastqc.zip"),
    log:
        str(LOG_DIR / "fastqc" / "{sample}_pretrim.log"),
    shell:
        """
        mkdir -p {QC_DIR}/pre_trim
        fastqc --outdir {QC_DIR}/pre_trim --threads {threads} --extract {input.fq} > {log} 2>&1
        # Rename if needed (fastqc appends _fastqc suffix to input stem)
        STEM=$(basename {input.fq} .fq.gz | sed 's/\.merged$//')
        for ext in html zip; do
            SRC="{QC_DIR}/pre_trim/${{STEM}}_fastqc.$ext"
            DST="{QC_DIR}/pre_trim/{wildcards.sample}_fastqc.$ext"
            [ -f "$SRC" ] && [ "$SRC" != "$DST" ] && mv "$SRC" "$DST" || true
        done
        """

rule fastqc_posttrim:
    input:
        fq = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
    output:
        html = str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.html"),
        zip  = str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.zip"),
    log:
        str(LOG_DIR / "fastqc" / "{sample}_posttrim.log"),
    shell:
        """
        mkdir -p {QC_DIR}/post_trim
        fastqc --outdir {QC_DIR}/post_trim --threads {threads} --extract {input.fq} > {log} 2>&1
        """
