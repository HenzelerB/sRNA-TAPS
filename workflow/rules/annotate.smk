# =============================================================================
# rules/annotate.smk — Gene annotation of TAPS call TSVs
# =============================================================================
#
# Runs after taps_call. For each biotype × sample TSV, intersects called
# cytosine positions with Ensembl GTF gene features and adds:
#
#   gene_name     — e.g. "hsa-miR-21-5p", "RNA45S5", "SNHG14"
#   gene_id       — Ensembl stable ID, e.g. "ENSG00000284190"
#   gene_biotype  — GTF biotype string, e.g. "miRNA", "rRNA", "lncRNA"
#
# Chromosome positions are preserved unchanged for IGV compatibility.
# Sites with no overlapping gene feature receive "." in all three columns.
# =============================================================================

rule annotate_taps_calls:
    """
    Annotate TAPS call TSVs with gene names from Ensembl GTF.
    Adds gene_name, gene_id, gene_biotype columns.
    Chromosome/position columns retained for IGV.
    """
    input:
        tsv = str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
        gtf = config["reference"]["gtf"],
    output:
        tsv = str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps_annotated.tsv"),
    log:
        str(LOG_DIR / "annotate" / "{sample}_{biotype}.log"),
    resources:
        mem_mb        = 8000,
        runtime       = 30,
        slurm_account = config["slurm"]["account"],
        slurm_extra   = "",
    script:
        "../scripts/annotate_taps.py"
