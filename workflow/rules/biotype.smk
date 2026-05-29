# =============================================================================
# rules/biotype.smk — RNA biotype splitting
# =============================================================================

rule biotype_split:
    """
    Split BAM reads into per-biotype BAMs using Ensembl GTF annotation.

    Priority order (highest wins when a read overlaps multiple features):
        miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other

    Uses binary-search interval lookup (~10x faster than linear scan).
    Handles chr-prefix mismatches between BAM and GTF automatically.
    """
    input:
        bam = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
        gtf = config["reference"]["gtf"],
    output:
        expand(
            str(BIOTYPE_DIR / "{biotype}" / "{{sample}}_{biotype}.sorted.bam"),
            biotype=BIOTYPES,
        ),
        summary = str(BIOTYPE_DIR / "{sample}_biotype_summary.txt"),
    params:
        out_dir = str(BIOTYPE_DIR),
        script  = str(Path(workflow.basedir).parent / "srnataps" / "annotate.py"),
    log:
        str(LOG_DIR / "biotype" / "{sample}.log"),
    threads: 4
    resources:
        mem_mb   = 30000,
        runtime  = 240,
    shell:
        """
        python {params.script} \
            --bam     {input.bam} \
            --gtf     {input.gtf} \
            --out_dir {params.out_dir} \
            --sample  {wildcards.sample} \
            > {log} 2>&1
        """


rule biotype_summarise:
    """Aggregate per-sample biotype summaries into a cross-sample table."""
    input:
        expand(str(BIOTYPE_DIR / "{sample}_biotype_summary.txt"), sample=SAMPLE_NAMES),
    output:
        str(BIOTYPE_DIR / "biotype_composition_all_samples.tsv"),
    params:
        in_dir = str(BIOTYPE_DIR),
        script = str(Path(workflow.basedir).parent / "srnataps" / "annotate.py"),
    log:
        str(LOG_DIR / "biotype" / "summarise.log"),
    shell:
        """
        python {params.script} summarise --in-dir {params.in_dir} > {log} 2>&1
        """
