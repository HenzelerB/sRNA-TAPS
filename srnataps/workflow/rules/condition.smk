# Condition-level pooling, dual-control testing, and evidence tiers.

rule pooled_condition_call:
    input:
        bams=condition_biotype_bams,
        fasta=config["reference"]["genome_fa"],
        snp_bed=lambda wc: str(SNP_DIR / f"sample_snps_{wc.cell_line}.bed"),
    output:
        tsv=str(POOLED_DIR / "{biotype}" / "{condition}_{cell_line}_pooled_{biotype}_taps.tsv"),
    params:
        script=str(SRNATAPS_SCRIPTS / "caller.py"),
        min_cov=CONDITION_CONFIG.get("pooled_min_coverage", 1),
        min_qual=config["calling"]["min_base_quality"],
        min_mapq=config["calling"]["min_mapping_quality"],
        bg_rate=config["calling"]["background_rate"],
        het_thresh=config["snp"]["het_threshold"],
        dbsnp=config["reference"].get("dbsnp_vcf", ""),
    threads: config.get("alignment", {}).get("threads", 8)
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        str(LOG_DIR / "condition" / "pooled_{condition}_{cell_line}_{biotype}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        DBSNP_ARG=""
        if [ -n "{params.dbsnp}" ] && [ -f "{params.dbsnp}" ]; then
            DBSNP_ARG="--dbsnp-vcf {params.dbsnp}"
        fi
        python {params.script} \
            --bam {input.bams} \
            --fasta {input.fasta} \
            --out {output.tsv} \
            --min-cov {params.min_cov} \
            --min-qual {params.min_qual} \
            --min-mapq {params.min_mapq} \
            --background-rate {params.bg_rate} \
            --sample-snp-bed {input.snp_bed} \
            --het-threshold {params.het_thresh} \
            --cell-line {wildcards.cell_line} \
            --threads {threads} \
            $DBSNP_ARG > {log} 2>&1
        """


rule pooled_control_contrast:
    input:
        treat=str(POOLED_DIR / "{biotype}" / "treat_{cell_line}_pooled_{biotype}_taps.tsv"),
        pb=str(POOLED_DIR / "{biotype}" / "pb_ctrl_{cell_line}_pooled_{biotype}_taps.tsv"),
        untreated=str(POOLED_DIR / "{biotype}" / "no_treat_{cell_line}_pooled_{biotype}_taps.tsv"),
    output:
        passing=str(CONTRAST_DIR / "{biotype}" / "treat_{cell_line}_contrast_{biotype}_taps.tsv"),
        audit=str(EVIDENCE_AUDIT_DIR / "pooled" / "{cell_line}_{biotype}_all_tested.tsv"),
    params:
        min_cov=CONDITION_CONFIG.get("test_min_coverage", 5),
        min_delta=CONDITION_CONFIG.get("minimum_delta", 0.1),
        max_padj=CONDITION_CONFIG.get("discovery_max_padj", 0.05),
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        str(LOG_DIR / "condition" / "contrast_{cell_line}_{biotype}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.passing}) $(dirname {output.audit}) $(dirname {log})
        python -m srnataps.contrast \
            --treat {input.treat} \
            --pb-control {input.pb} \
            --no-treat {input.untreated} \
            --out {output.passing} \
            --all-out {output.audit} \
            --min-treat-coverage {params.min_cov} \
            --min-control-coverage {params.min_cov} \
            --min-delta {params.min_delta} \
            --max-padj {params.max_padj} > {log} 2>&1
        """


rule replicate_condition_call:
    input:
        pooled_audit=str(EVIDENCE_AUDIT_DIR / "pooled" / "{cell_line}_{biotype}_all_tested.tsv"),
        bams=all_biotype_bams,
        samples="samples.tsv",
    output:
        passing=str(REPLICATE_DIR / "{biotype}" / "treat_{cell_line}_replicate_{biotype}_taps.tsv"),
        audit=str(EVIDENCE_AUDIT_DIR / "replicate" / "{cell_line}_{biotype}_all_tested.tsv"),
    params:
        min_reps=CONDITION_CONFIG.get("minimum_replicates", 3),
        min_delta=CONDITION_CONFIG.get("minimum_delta", 0.1),
        max_padj=CONDITION_CONFIG.get("discovery_max_padj", 0.05),
        min_cov=CONDITION_CONFIG.get("test_min_coverage", 5),
        min_qual=config["calling"]["min_base_quality"],
        min_mapq=config["calling"]["min_mapping_quality"],
    resources:
        mem_mb=12000,
        runtime=480,
    log:
        str(LOG_DIR / "condition" / "replicate_{cell_line}_{biotype}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.passing}) $(dirname {output.audit}) $(dirname {log})
        python -m srnataps.replicate_contrast \
            --audit {input.pooled_audit} \
            --calls-dir {CALLS_DIR} \
            --bam-dir {BIOTYPE_DIR} \
            --samples-tsv {input.samples} \
            --biotype {wildcards.biotype} \
            --cell-line {wildcards.cell_line} \
            --out {output.passing} \
            --all-out {output.audit} \
            --statistical-test beta_binomial \
            --minimum-replicates {params.min_reps} \
            --minimum-delta {params.min_delta} \
            --maximum-padj {params.max_padj} \
            --minimum-control-coverage {params.min_cov} \
            --minimum-base-quality {params.min_qual} \
            --minimum-mapping-quality {params.min_mapq} \
            --minimum-sample-coverage 1 > {log} 2>&1
        """


rule stringent_condition_call:
    input:
        audit=str(EVIDENCE_AUDIT_DIR / "replicate" / "{cell_line}_{biotype}_all_tested.tsv"),
    output:
        tsv=str(STRINGENT_DIR / "{biotype}" / "treat_{cell_line}_stringent_{biotype}_taps.tsv"),
    params:
        max_padj=CONDITION_CONFIG.get("stringent_max_padj", 1e-20),
        min_delta=CONDITION_CONFIG.get("minimum_delta", 0.1),
        min_cov=CONDITION_CONFIG.get("stringent_min_coverage", 5),
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        str(LOG_DIR / "condition" / "stringent_{cell_line}_{biotype}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -m srnataps.tier \
            --audit {input.audit} \
            --out {output.tsv} \
            --maximum-padj {params.max_padj} \
            --minimum-delta {params.min_delta} \
            --minimum-coverage {params.min_cov} > {log} 2>&1
        """


rule evidence_tier_summary:
    input:
        pooled=expand(
            str(POOLED_DIR / "{biotype}" / "treat_{cell_line}_pooled_{biotype}_taps.tsv"),
            cell_line=CELL_LINES, biotype=BIOTYPES,
        ),
        contrast=expand(
            str(CONTRAST_DIR / "{biotype}" / "treat_{cell_line}_contrast_{biotype}_taps.tsv"),
            cell_line=CELL_LINES, biotype=BIOTYPES,
        ),
        replicate=expand(
            str(REPLICATE_DIR / "{biotype}" / "treat_{cell_line}_replicate_{biotype}_taps.tsv"),
            cell_line=CELL_LINES, biotype=BIOTYPES,
        ),
        stringent=expand(
            str(STRINGENT_DIR / "{biotype}" / "treat_{cell_line}_stringent_{biotype}_taps.tsv"),
            cell_line=CELL_LINES, biotype=BIOTYPES,
        ),
    output:
        tsv=str(REPORT_DIR / "tables" / "evidence_tier_summary.tsv"),
    resources:
        mem_mb=2000,
        runtime=30,
    log:
        str(LOG_DIR / "condition" / "evidence_tier_summary.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -m srnataps.evidence_summary \
            --cell-lines {CELL_LINES} \
            --biotypes {BIOTYPES} \
            --pooled-dir {POOLED_DIR} \
            --contrast-dir {CONTRAST_DIR} \
            --replicate-dir {REPLICATE_DIR} \
            --stringent-dir {STRINGENT_DIR} \
            --out {output.tsv} > {log} 2>&1
        """
