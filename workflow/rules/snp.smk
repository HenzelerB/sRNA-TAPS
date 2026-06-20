# =============================================================================
# rules/snp.smk — SNP blacklist construction (per cell line)
# =============================================================================

rule merge_notreat_bams:
    """Merge no-treat replicate BAMs per cell line for SNP calling."""
    input:
        bams = get_notreat_bams,
    output:
        bam = str(SNP_DIR / "notreat_{cell_line}_merged.bam"),
        bai = str(SNP_DIR / "notreat_{cell_line}_merged.bam.bai"),
    log:
        str(LOG_DIR / "snp" / "merge_notreat_{cell_line}.log"),
    shell:
        """
        mkdir -p {SNP_DIR}
        if [ $(echo {input.bams} | wc -w) -gt 1 ]; then
            samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        else
            cp {input.bams} {output.bam}
        fi
        samtools index {output.bam} {output.bai}
        """

rule snp_blacklist:
    """
    Build per-cell-line SNP blacklist from no-treat BAMs.
    C→T / G→A at AF >= min_af in no-treat = germline variant, not modification.
    """
    input:
        bam   = str(SNP_DIR / "notreat_{cell_line}_merged.bam"),
        fasta = config["reference"]["genome_fa"],
    output:
        bed = str(SNP_DIR / "sample_snps_{cell_line}.bed"),
    params:
        script  = str(SRNATAPS_SCRIPTS / "snp.py"),
        min_af  = config["snp"]["min_af"],
        min_cov = config["snp"]["min_cov"],
    log:
        str(LOG_DIR / "snp" / "blacklist_{cell_line}.log"),
    shell:
        """
        python {params.script} \
            --bam       {input.bam} \
            --fasta     {input.fasta} \
            --out       {output.bed} \
            --min-af    {params.min_af} \
            --min-cov   {params.min_cov} \
            --cell-line {wildcards.cell_line} \
            > {log} 2>&1
        """
