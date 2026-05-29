# =============================================================================
# rules/snp.smk — SNP blacklist construction (three-layer filtering)
# =============================================================================

def get_notreat_bams(wc):
    """Return all no-treat BAMs for a given cell line."""
    notreat = SAMPLES[
        (SAMPLES["condition"] == "no_treat") &
        (SAMPLES["cell_line"] == wc.cell_line)
    ]["sample"].tolist()
    return [str(ALIGN_DIR / f"{s}.sorted.bam") for s in notreat]


rule merge_notreat_bams:
    """
    Merge no-treat replicate BAMs per cell line for SNP calling.
    Higher coverage from merged replicates = more confident SNP calls.
    """
    input:
        bams = get_notreat_bams,
    output:
        bam = str(SNP_DIR / "notreat_{cell_line}_merged.bam"),
        bai = str(SNP_DIR / "notreat_{cell_line}_merged.bam.bai"),
    log:
        str(LOG_DIR / "snp" / "merge_notreat_{cell_line}.log"),
    threads: 4
    resources:
        mem_mb   = 16000,
        runtime  = 120,
    shell:
        """
        mkdir -p {SNP_DIR}
        if [ $(echo {input.bams} | wc -w) -gt 1 ]; then
            samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        else
            cp {input.bams} {output.bam}
        fi
        samtools index {output.bam}
        """


rule build_snp_blacklist:
    """
    Build per-cell-line SNP blacklist from no-treat BAMs.

    The no-treat condition has no TAPS chemistry applied.
    Any C→T or G→A at AF >= min_af is a germline variant, not modification.

    HEK293 and Caco2 are immortalised lines with distinct mutation profiles
    not fully captured in population dbSNP — this per-sample approach catches
    cell-line-specific variants that dbSNP alone would miss.

    Output BED columns: chrom, start, end, ref, alt, allele_freq, coverage, cell_line
    """
    input:
        bam   = str(SNP_DIR / "notreat_{cell_line}_merged.bam"),
        fasta = config["reference"]["genome_fa"],
    output:
        bed = str(SNP_DIR / "sample_snps_{cell_line}.bed"),
    params:
        script  = str(Path(workflow.basedir).parent / "srnataps" / "snp.py"),
        min_af  = config["snp"]["min_af"],
        min_cov = config["snp"]["min_cov"],
    log:
        str(LOG_DIR / "snp" / "blacklist_{cell_line}.log"),
    threads: 1
    resources:
        mem_mb   = 30000,
        runtime  = 240,
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
