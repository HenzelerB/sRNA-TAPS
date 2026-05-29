# =============================================================================
# rules/call.smk — TAPS m5C modification calling
# =============================================================================

def get_snp_bed(wc):
    """Return the SNP blacklist BED for the cell line of this sample."""
    cell_line = SAMPLES.loc[wc.sample, "cell_line"]
    return str(SNP_DIR / f"sample_snps_{cell_line}.bed")


def get_min_cov(wc):
    """Return per-biotype minimum coverage threshold."""
    return config["calling"]["min_coverage"].get(wc.biotype, 5)


rule taps_call:
    """
    TAPS m5C modification calling per biotype per sample.

    Chemistry:
        TAPS: m5C and 5hmC → T (TET oxidation + pyridine borane reduction)
        Unmodified C stays as C.
        C→T in reads at a reference-C = modification (opposite of bisulfite).
        Reverse strand: G→A = modification at minus-strand C.

    SNP filtering (three layers, applied BEFORE counting):
        Layer 1: dbSNP common C→T/G→A variants (AF >= 0.01)
        Layer 2: Cell-line-specific SNPs from no-treat BAMs
        Layer 3: Heterozygosity flag (no-treat C→T rate >= het_threshold)
        SNP sites are EXCLUDED before pileup counting, ensuring the
        BH FDR correction pool contains only genuine modification candidates.

    Multi-mapper weighting:
        Bowtie1 writes XA:i:N (number of alignments).
        Reads mapping to N loci each contribute weight = 1/N.
        Essential for tRNA (600+ near-identical copies in hg38).

    Output columns:
        chrom, start, end, context, mod_count, unmod_count, coverage,
        mod_rate, pvalue, padj, snp_flag
    """
    input:
        bam     = str(BIOTYPE_DIR / "{biotype}" / "{sample}_{biotype}.sorted.bam"),
        fasta   = config["reference"]["genome_fa"],
        snp_bed = get_snp_bed,
        snp_idx = lambda wc: str(SNP_DIR / f"sample_snps_{SAMPLES.loc[wc.sample, 'cell_line']}.bed"),
    output:
        tsv = str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
    params:
        script      = str(Path(workflow.basedir).parent / "srnataps" / "caller.py"),
        min_cov     = get_min_cov,
        min_qual    = config["calling"]["min_base_quality"],
        min_mapq    = config["calling"]["min_mapping_quality"],
        bg_rate     = config["calling"]["background_rate"],
        het_thresh  = config["snp"]["het_threshold"],
        dbsnp       = config["reference"].get("dbsnp_vcf", ""),
        cell_line   = lambda wc: SAMPLES.loc[wc.sample, "cell_line"],
    log:
        str(LOG_DIR / "call" / "{sample}_{biotype}.log"),
    threads: 4
    resources:
        mem_mb   = 40000,
        runtime  = 480,
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
