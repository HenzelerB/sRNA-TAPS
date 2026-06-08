# =============================================================================
# rules/qc.smk — FastQC and Trim Galore
# =============================================================================

rule fastqc_pretrim:
    """FastQC on raw merged FASTQs before adapter trimming."""
    input:
        fq = lambda wc: SAMPLES.loc[wc.sample, "fastq"],
    output:
        html = str(QC_DIR / "pre_trim" / "{sample}_fastqc.html"),
        zip  = str(QC_DIR / "pre_trim" / "{sample}_fastqc.zip"),
    params:
        outdir = str(QC_DIR / "pre_trim"),
    log:
        str(LOG_DIR / "fastqc" / "pretrim_{sample}.log"),
    threads: 4
    resources:
        mem_mb   = 8000,
        runtime  = 60,
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --outdir {params.outdir} --threads {threads} --extract {input.fq} \
            > {log} 2>&1
        # Rename output to match expected filenames (strip .merged suffix)
        RAW=$(basename {input.fq})
        RAW="${{RAW%.fastq.gz}}"
        RAW="${{RAW%.fq.gz}}"
        SRC_HTML={params.outdir}/${{RAW}}_fastqc.html
        SRC_ZIP={params.outdir}/${{RAW}}_fastqc.zip
        if [ -f "$SRC_HTML" ] && [ "$SRC_HTML" != "{output.html}" ]; then
            mv "$SRC_HTML" {output.html}
            mv "$SRC_ZIP"  {output.zip}
        fi
        """


rule trim_galore:
    """
    Adapter trimming with Trim Galore using TruSeq small RNA adapter.

    TAPS note: do NOT use --rrbs or any bisulfite mode.
    TAPS preserves unmodified C — bisulfite modes would misinterpret chemistry.
    --small_rna activates the TruSeq SR adapter (TGGAATTCTCGGGTGCCAAGG)
    and sets --length 18 (discard reads < 18 nt after trimming).
    """
    input:
        fq = lambda wc: SAMPLES.loc[wc.sample, "fastq"],
    output:
        fq      = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        report  = str(TRIM_DIR / "{sample}_trimming_report.txt"),
    params:
        outdir    = str(TRIM_DIR),
        adapter   = config["trimming"]["adapter"],
        min_len   = config["trimming"]["min_length"],
        max_len   = config["trimming"]["max_length"],
    log:
        str(LOG_DIR / "trim" / "{sample}.log"),
    threads: 4
    resources:
        mem_mb   = 8000,
        runtime  = 120,
    shell:
        """
        trim_galore \
            --small_rna \
            --cores {threads} \
            -o {params.outdir} \
            {input.fq} \
            > {log} 2>&1

        # Rename to consistent name (use exact match to avoid greedy glob)
        INPUT_BASE=$(basename {input.fq})
        INPUT_STEM="${{INPUT_BASE%.fastq.gz}}"
        INPUT_STEM="${{INPUT_STEM%.fq.gz}}"
        TRIMMED=$(ls {params.outdir}/${{INPUT_STEM}}_trimmed.fq.gz 2>/dev/null || true)
        REPORT=$(ls  {params.outdir}/${{INPUT_STEM}}_trimming_report.txt 2>/dev/null || true)
        if [[ -n "$TRIMMED" && "$TRIMMED" != "{output.fq}" ]]; then
            mv "$TRIMMED" "{output.fq}"
        fi
        if [[ -n "$REPORT" && "$REPORT" != "{output.report}" ]]; then
            mv "$REPORT" "{output.report}"
        fi
        """


rule fastqc_posttrim:
    """FastQC on trimmed FASTQs."""
    input:
        fq = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
    output:
        html = str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.html"),
        zip  = str(QC_DIR / "post_trim" / "{sample}_trimmed_fastqc.zip"),
    params:
        outdir = str(QC_DIR / "post_trim"),
    log:
        str(LOG_DIR / "fastqc" / "posttrim_{sample}.log"),
    threads: 4
    resources:
        mem_mb   = 8000,
        runtime  = 60,
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --outdir {params.outdir} --threads {threads} --extract {input.fq} \
            > {log} 2>&1
        # Rename output to match expected filenames (strip .merged suffix)
        RAW=$(basename {input.fq})
        RAW="${{RAW%.fastq.gz}}"
        RAW="${{RAW%.fq.gz}}"
        SRC_HTML={params.outdir}/${{RAW}}_fastqc.html
        SRC_ZIP={params.outdir}/${{RAW}}_fastqc.zip
        if [ -f "$SRC_HTML" ] && [ "$SRC_HTML" != "{output.html}" ]; then
            mv "$SRC_HTML" {output.html}
            mv "$SRC_ZIP"  {output.zip}
        fi
        """
