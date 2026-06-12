# =============================================================================
# rules/benchmark.smk — rastair, asTair, Bismark benchmarking
# =============================================================================

rule rastair_cpg:
    """rastair 2.1.1 — CpG-only calling on whole-genome Bowtie1 BAMs."""
    input:
        bam  = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai  = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
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
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -f "{input.fai}" ]; then
            samtools faidx {input.fasta}
        fi
        rastair call \
            -r            {input.fasta} \
            --bed         {output.bed} \
            --unpaired \
            -@            {threads} \
            --v-min-depth {params.min_depth} \
            -Q            {params.min_baseq} \
            -q            {params.min_mapq} \
            --cpgs-only \
            --no-ml \
            {input.bam} \
            > {log} 2>&1
        """

rule rastair_all:
    """rastair 2.1.1 — all-context calling on whole-genome Bowtie1 BAMs."""
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
    shell:
        """
        mkdir -p {params.outdir}
        rastair call \
            -r            {input.fasta} \
            --bed         {output.bed} \
            --unpaired \
            -@            {threads} \
            --v-min-depth {params.min_depth} \
            -Q            {params.min_baseq} \
            -q            {params.min_mapq} \
            --no-ml \
            {input.bam} \
            > {log} 2>&1
        """

rule astair_call:
    """
    asTair 3.3.2 — all-context calling on biotype-split Bowtie1 BAMs.
    -m mCtoT: TAPS chemistry
    -co all : all C contexts (essential for non-CpG m5C in tRNA/rRNA)
    Removes existing .mods.gz before running to avoid 'file exists' error.
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
    resources:
        mem_mb  = lambda wc, input: est_mem(
            10000, input.bam,
            scale=15, floor_mb=12000, ceil_mb=80000
        ),
        runtime = 720,
    log:
        str(LOG_DIR / "benchmark" / "astair_{sample}_{biotype}.log"),
    shell:
        """
        mkdir -p {params.outdir}

        # Remove any existing mods files to avoid asTair "file exists" error
        rm -f {params.outdir}/$(basename {input.bam} .bam)_mCtoT_{params.context}.mods.gz

        READ_COUNT=$(samtools view -c -F 4 {input.bam})
        if [ "$READ_COUNT" -lt 50 ]; then
            echo "[SKIP] Too few reads ($READ_COUNT) for {wildcards.sample} {wildcards.biotype}"
            touch {output.done}
            exit 0
        fi

        export PYTHONNOUSERSITE=1
        /opt/apps/conda/bhenzeler/envs/astair_env/bin/astair call \
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

rule bismark_index:
    """Build Bismark bowtie2 genome index — run once. Uses TMPDIR in index dir to avoid /tmp overflow."""
    input:
        fa = config["reference"]["genome_fa"],
    output:
        done = str(BENCH_DIR / "bismark_index" / ".index_complete"),
    params:
        index_dir = str(BENCH_DIR / "bismark_index"),
    log:
        str(LOG_DIR / "benchmark" / "bismark_index.log"),
    shell:
        """
        mkdir -p {params.index_dir}
        cp {input.fa} {params.index_dir}/
        export TMPDIR={params.index_dir}
        bismark_genome_preparation \
            --bowtie2 \
            --large-index \
            --parallel 2 \
            {params.index_dir} \
            > {log} 2>&1
        touch {output.done}
        """

rule bismark_align:
    """Bismark alignment + methylation extraction per sample."""
    input:
        fq   = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        done = str(BENCH_DIR / "bismark_index" / ".index_complete"),
        fasta = config["reference"]["genome_fa"],
    output:
        done = str(BENCH_DIR / "bismark" / "{sample}" / ".done"),
    params:
        outdir    = str(BENCH_DIR / "bismark" / "{sample}"),
        index_dir = str(BENCH_DIR / "bismark_index"),
        min_mapq  = config["benchmark"]["tools"]["bismark"]["min_mapq"],
        min_baseq = config["benchmark"]["tools"]["bismark"]["min_baseq"],
    resources:
        mem_mb  = lambda wc, input: est_mem(
            14000, input.fq,
            scale=15, floor_mb=16000, ceil_mb=60000
        ),
        runtime = 480,
    log:
        align   = str(LOG_DIR / "benchmark" / "bismark_align_{sample}.log"),
        extract = str(LOG_DIR / "benchmark" / "bismark_extract_{sample}.log"),
    shell:
        """
        mkdir -p {params.outdir}

        bismark \
            --bowtie2 \
            --non_directional \
            --score_min L,-1,-0.6 \
            --rdg 100,100 \
            --rfg 100,100 \
            -p {threads} \
            --bam \
            -o {params.outdir} \
            {params.index_dir} \
            {input.fq} \
            > {log.align} 2>&1

        BAM=$(ls {params.outdir}/*.bam 2>/dev/null | grep -v bismark_bt2 | head -1 || true)
        if [ -z "$BAM" ]; then
            BAM=$(ls {params.outdir}/*.bam 2>/dev/null | head -1 || true)
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
