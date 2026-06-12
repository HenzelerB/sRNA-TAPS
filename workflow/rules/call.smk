# =============================================================================
# rules/call.smk — TAPS m5C modification calling per biotype per sample
# =============================================================================

rule taps_call:
    """
    TAPS m5C calling per biotype per sample.
    Three-layer SNP filtering: dbSNP + cell-line-specific + heterozygosity.
    Multi-mapper weighting via XA:i:N tag (essential for tRNA).
    Skips biotypes with < 50 reads (writes empty TSV).
    """
    input:
        bam     = str(BIOTYPE_DIR / "{biotype}" / "{sample}_{biotype}.sorted.bam"),
        fasta   = config["reference"]["genome_fa"],
        snp_bed = get_snp_bed,
    output:
        tsv = str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
    params:
        script    = str(SRNATAPS_SCRIPTS / "caller.py"),
        min_cov   = get_min_cov,
        min_qual  = config["calling"]["min_base_quality"],
        min_mapq  = config["calling"]["min_mapping_quality"],
        bg_rate   = config["calling"]["background_rate"],
        het_thresh = config["snp"]["het_threshold"],
        cell_line = lambda wc: get_cell_line(wc.sample),
        dbsnp     = config["reference"].get("dbsnp_vcf", ""),
    resources:
        mem_mb  = lambda wc, input: est_mem(
            6000, input.bam,
            scale=12, floor_mb=8000, ceil_mb=64000
        ),
        runtime = 480,
    log:
        str(LOG_DIR / "call" / "{sample}_{biotype}.log"),
    shell:
        """
        mkdir -p {CALLS_DIR}/{wildcards.biotype}

        # Skip if BAM has too few reads
        READ_COUNT=$(samtools view -c -F 4 {input.bam})
        if [ "$READ_COUNT" -lt 50 ]; then
            echo "[SKIP] Too few reads ($READ_COUNT) for {wildcards.sample} {wildcards.biotype}"
            touch {output.tsv}
            exit 0
        fi

        DBSNP_ARG=""
        if [ -n "{params.dbsnp}" ] && [ -f "{params.dbsnp}" ]; then
            DBSNP_ARG="--dbsnp-vcf {params.dbsnp}"
        fi

        python {params.script} \
            --bam               {input.bam} \
            --fasta             {input.fasta} \
            --out               {output.tsv} \
            --min-cov           {params.min_cov} \
            --min-qual          {params.min_qual} \
            --min-mapq          {params.min_mapq} \
            --background-rate   {params.bg_rate} \
            --sample-snp-bed    {input.snp_bed} \
            --het-threshold     {params.het_thresh} \
            --cell-line         {params.cell_line} \
            --threads           {threads} \
            $DBSNP_ARG \
            > {log} 2>&1
        """
