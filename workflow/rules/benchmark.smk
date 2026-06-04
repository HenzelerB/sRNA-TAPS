# =============================================================================
# rules/benchmark.smk — Benchmarking: rastair, asTair, Bismark
# Only runs when --benchmark flag is set or config benchmark.enabled = true
# =============================================================================

# ── rastair ───────────────────────────────────────────────────────────────────
rule rastair_call_cpg:
    """
    rastair 2.1.1 CpG-only methylation calling on Bowtie1-aligned BAMs.
    Output: 08.benchmark/rastair/cpg/<sample>/
    Used for CpG-specific concordance analysis vs sRNA-TAPS.
    """
    input:
        bam   = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai   = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
        fasta = config["reference"]["genome_fa"],
        fai   = config["reference"]["genome_fa"] + ".fai",
    output:
        bed = str(BENCH_DIR / "rastair" / "cpg" / "{sample}" / "{sample}.bed.gz"),
    params:
        outdir    = str(BENCH_DIR / "rastair" / "cpg" / "{sample}"),
        min_depth = config["benchmark"]["tools"]["rastair"]["min_depth"],
        min_mapq  = config["benchmark"]["tools"]["rastair"]["min_mapq"],
        min_baseq = config["benchmark"]["tools"]["rastair"]["min_baseq"],
    log:
        str(LOG_DIR / "benchmark" / "rastair_cpg_{sample}.log"),
    threads: 8
    resources:
        mem_mb   = 32000,
        runtime  = 240,
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -f "{input.fai}" ]; then
            samtools faidx {input.fasta}
        fi
        rastair call \
            -r              {input.fasta} \
            --bed           {output.bed} \
            --unpaired \
            -@              {threads} \
            --v-min-depth   {params.min_depth} \
            -Q              {params.min_baseq} \
            -q              {params.min_mapq} \
            --cpgs-only \
            --no-ml \
            {input.bam} \
            > {log} 2>&1
        """


rule rastair_call_all:
    """
    rastair 2.1.1 all-context methylation calling on Bowtie1-aligned BAMs.
    Output: 08.benchmark/rastair/all/<sample>/
    Used for full concordance analysis vs sRNA-TAPS (all cytosine contexts).
    Note: rastair reports all positions it can call; non-CpG sites included.
    """
    input:
        bam   = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai   = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
        fasta = config["reference"]["genome_fa"],
        fai   = config["reference"]["genome_fa"] + ".fai",
    output:
        bed = str(BENCH_DIR / "rastair" / "all" / "{sample}" / "{sample}.bed.gz"),
    params:
        outdir    = str(BENCH_DIR / "rastair" / "all" / "{sample}"),
        min_depth = config["benchmark"]["tools"]["rastair"]["min_depth"],
        min_mapq  = config["benchmark"]["tools"]["rastair"]["min_mapq"],
        min_baseq = config["benchmark"]["tools"]["rastair"]["min_baseq"],
    log:
        str(LOG_DIR / "benchmark" / "rastair_all_{sample}.log"),
    threads: 8
    resources:
        mem_mb   = 32000,
        runtime  = 240,
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -f "{input.fai}" ]; then
            samtools faidx {input.fasta}
        fi
        rastair call \
            -r              {input.fasta} \
            --bed           {output.bed} \
            --unpaired \
            -@              {threads} \
            --v-min-depth   {params.min_depth} \
            -Q              {params.min_baseq} \
            -q              {params.min_mapq} \
            --no-ml \
            {input.bam} \
            > {log} 2>&1
        """


# ── asTair ────────────────────────────────────────────────────────────────────
rule astair_call:
    """
    asTair 3.3.2 modification calling on Bowtie1 biotype-split BAMs.
    -m mCtoT : TAPS chemistry (modified C → T)
    -co all  : all C contexts (CpG + CHH + CHG) — essential for non-CpG m5C in tRNA/rRNA
    -se      : single-end reads
    Fully matched comparison: same Bowtie1 BAMs, same biotype separation.
    """
    input:
        bam   = str(BIOTYPE_DIR / "{biotype}" / "{sample}_{biotype}.sorted.bam"),
        fasta = config["reference"]["genome_fa"],
    output:
        done = str(BENCH_DIR / "astair" / "{biotype}" / ".done_{sample}"),
    params:
        outdir    = str(BENCH_DIR / "astair" / "{biotype}"),
        min_depth = config["benchmark"]["tools"]["astair"]["min_depth"],
        min_mapq  = config["benchmark"]["tools"]["astair"]["min_mapq"],
        min_baseq = config["benchmark"]["tools"]["astair"]["min_baseq"],
        context   = config["benchmark"]["tools"]["astair"]["context"],
    log:
        str(LOG_DIR / "benchmark" / "astair_{sample}_{biotype}.log"),
    threads: 8
    resources:
        mem_mb   = 32000,
        runtime  = 240,
    shell:
        """
        mkdir -p {params.outdir}

        READ_COUNT=$(samtools view -c -F 4 {input.bam})
        if [ "$READ_COUNT" -lt 50 ]; then
            echo "[SKIP] Too few reads for {wildcards.sample} {wildcards.biotype}"
            touch {output.done}
            exit 0
        fi

        export PYTHONNOUSERSITE=1
        astair call \
            -i  {input.bam} \
            -f  {input.fasta} \
            -m  mCtoT \
            -co {params.context} \
            -se \
            -bq {params.min_baseq} \
            -mq {params.min_mapq} \
            -t  {threads} \
            -d  {params.outdir} \
            -z \
            > {log} 2>&1

        touch {output.done}
        """


# ── Bismark ───────────────────────────────────────────────────────────────────
rule bismark_index:
    """
    Build Bismark bisulfite genome index (bowtie2, once per project).
    Note: Bismark v0.24.2 does not support bowtie1 index build.
    Bismark is used as a chemistry contrast reference — its C→T interpretation
    is the INVERSE of TAPS. In 09_compare.py: taps_equiv = 1 - bismark_rate.
    """
    input:
        fa = config["reference"]["genome_fa"],
    output:
        done = str(BENCH_DIR / "bismark_index" / ".index_complete"),
    params:
        index_dir = str(BENCH_DIR / "bismark_index"),
    log:
        str(LOG_DIR / "benchmark" / "bismark_index.log"),
    threads: 8
    resources:
        mem_mb   = 64000,
        runtime  = 480,
    shell:
        """
        mkdir -p {params.index_dir}
        cp {input.fa} {params.index_dir}/
        bismark_genome_preparation \
            --bowtie2 \
            --large-index \
            --parallel 4 \
            {params.index_dir} \
            > {log} 2>&1
        touch {output.done}
        """


rule bismark_align_extract:
    """
    Bismark alignment + methylation extraction per sample.
    Chemistry inversion: bismark_rate = unmethylated fraction in bisulfite.
    09_compare.py applies (1 - bismark_rate) before comparison with TAPS calls.
    """
    input:
        fq    = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        index_done = str(BENCH_DIR / "bismark_index" / ".index_complete"),
    output:
        done = str(BENCH_DIR / "bismark" / "{sample}" / ".done"),
    params:
        outdir    = str(BENCH_DIR / "bismark" / "{sample}"),
        index_dir = str(BENCH_DIR / "bismark_index"),
    log:
        align   = str(LOG_DIR / "benchmark" / "bismark_align_{sample}.log"),
        extract = str(LOG_DIR / "benchmark" / "bismark_extract_{sample}.log"),
    threads: 8
    resources:
        mem_mb   = 48000,
        runtime  = 480,
    shell:
        """
        mkdir -p {params.outdir}

        bismark \
            --bowtie2 \
            --non_directional \
            --score_min L,-1,-0.6 \
            --rdg 100,100 \
            --rfg 100,100 \
            --output_dir "{params.outdir}" \
            "{params.index_dir}" \
            "{input.fq}" \
            > {log.align} 2>&1

        BAM=$(ls {params.outdir}/*.bam 2>/dev/null | head -1)
        if [ -z "$BAM" ]; then
            echo "ERROR: No BAM produced for {wildcards.sample}"
            exit 1
        fi

        samtools sort -@ 4 -o {params.outdir}/{wildcards.sample}.sorted.bam "$BAM"
        samtools index {params.outdir}/{wildcards.sample}.sorted.bam
        rm -f "$BAM"

        bismark_methylation_extractor \
            --single-end \
            --CX_context \
            --comprehensive \
            --bedGraph \
            -o {params.outdir} \
            {params.outdir}/{wildcards.sample}.sorted.bam \
            > {log.extract} 2>&1

        touch {output.done}
        """
