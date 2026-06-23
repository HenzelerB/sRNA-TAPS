# =============================================================================
# rules/biotype.smk — Per-biotype BAM splitting
# =============================================================================

SRNATAPS_SCRIPTS = Path(workflow.basedir).parent

rule biotype_split:
    """
    Annotate and split BAM by RNA biotype using Ensembl GTF.
    Priority: miRNA > tRNA > piRNA > snoRNA > snRNA > rRNA > lncRNA > other
    Creates empty sorted BAMs for biotypes with 0 reads (piRNA in most samples).
    """
    input:
        bam = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
        gtf = config["reference"]["gtf"],
    output:
        expand(
            str(BIOTYPE_DIR / "{biotype}" / "{{sample}}_{biotype}.sorted.bam"),
            biotype=BIOTYPES
        ),
        summary = str(BIOTYPE_DIR / "{sample}_biotype_summary.txt"),
    log:
        str(LOG_DIR / "biotype" / "{sample}.log"),
    params:
        script  = str(SRNATAPS_SCRIPTS / "annotate.py"),
        out_dir = str(BIOTYPE_DIR),
    shell:
        """
        python {params.script} \
            --bam     {input.bam} \
            --gtf     {input.gtf} \
            --out_dir {params.out_dir} \
            --sample  {wildcards.sample} \
            > {log} 2>&1
        """
