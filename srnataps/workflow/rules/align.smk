# =============================================================================
# rules/align.smk - genome setup, Bowtie1 index, and small-RNA alignment
# =============================================================================

ENSEMBL_RELEASE = REFERENCE_DETAILS["ensembl_release"]
ENSEMBL_FA_URL = REFERENCE_DETAILS["genome_url"]
ENSEMBL_GTF_URL = REFERENCE_DETAILS["gtf_url"]

ALIGNMENT_STRATEGY = str(
    config.get("alignment", {}).get("strategy", "three_letter")
).lower()
if ALIGNMENT_STRATEGY not in {"three_letter", "standard"}:
    raise ValueError("alignment.strategy must be 'three_letter' or 'standard'")

THREE_LETTER_DIR = GENOME_DIR / "three_letter"
SRNATAPS_ROOT = Path(workflow.basedir).parent
REFERENCE_STEM = REFERENCE_DETAILS["index_stem"]
THREE_LETTER_C2T_FA = THREE_LETTER_DIR / f"{REFERENCE_STEM}_C2T.fa"
THREE_LETTER_G2A_FA = THREE_LETTER_DIR / f"{REFERENCE_STEM}_G2A.fa"
THREE_LETTER_C2T_INDEX = THREE_LETTER_DIR / f"{REFERENCE_STEM}_C2T"
THREE_LETTER_G2A_INDEX = THREE_LETTER_DIR / f"{REFERENCE_STEM}_G2A"


def alignment_index_markers(wildcards):
    if ALIGNMENT_STRATEGY == "three_letter":
        return [
            str(THREE_LETTER_DIR / ".c2t_index_complete"),
            str(THREE_LETTER_DIR / ".g2a_index_complete"),
        ]
    return [str(GENOME_DIR / ".bowtie1_index_complete")]


def bowtie_alignment_args():
    """Build the configured Bowtie1 small-RNA alignment policy."""
    alignment = config.get("alignment", {})
    mode = str(alignment.get("mode", "seed")).lower()
    mismatches = int(alignment.get("mismatches", 1))
    if not 0 <= mismatches <= 3:
        raise ValueError("alignment.mismatches must be between 0 and 3")

    if mode == "seed":
        seed_length = int(alignment.get("seed_length", 10))
        max_quality_sum = int(alignment.get("max_mismatch_quality", 100))
        return f"-n {mismatches} -l {seed_length} -e {max_quality_sum}"
    if mode in {"end_to_end", "v"}:
        return f"-v {mismatches}"
    raise ValueError("alignment.mode must be 'seed' or 'end_to_end'")


def bowtie_strand_arg():
    strand = str(config.get("alignment", {}).get("strand", "both")).lower()
    if strand == "both":
        return ""
    if strand in {"forward", "fw"}:
        return "--norc"
    if strand in {"reverse", "rc"}:
        return "--nofw"
    raise ValueError("alignment.strand must be 'both', 'forward', or 'reverse'")


rule download_genome_fa:
    output:
        fa=config["reference"]["genome_fa"],
    params:
        url=ENSEMBL_FA_URL,
        fa_gz=config["reference"]["genome_fa"] + ".gz",
        outdir=str(GENOME_DIR),
    log:
        str(LOG_DIR / "align" / "download_genome_fa.log"),
    shell:
        """
        mkdir -p {params.outdir}
        wget -q --show-progress -O {params.fa_gz} {params.url} > {log} 2>&1
        gunzip -f {params.fa_gz} >> {log} 2>&1
        """


rule download_gtf:
    output:
        gtf=config["reference"]["gtf"],
    params:
        url=ENSEMBL_GTF_URL,
        gtf_gz=config["reference"]["gtf"] + ".gz",
        outdir=str(GENOME_DIR),
    log:
        str(LOG_DIR / "align" / "download_gtf.log"),
    shell:
        """
        mkdir -p {params.outdir}
        wget -q --show-progress -O {params.gtf_gz} {params.url} > {log} 2>&1
        gunzip -f {params.gtf_gz} >> {log} 2>&1
        """


rule index_genome_fa:
    input:
        fa=config["reference"]["genome_fa"],
    output:
        fai=config["reference"]["genome_fa"] + ".fai",
    log:
        str(LOG_DIR / "align" / "index_genome_fa.log"),
    shell:
        """
        samtools faidx {input.fa} > {log} 2>&1
        """


rule bowtie1_index:
    input:
        fa=config["reference"]["genome_fa"],
    output:
        done=str(GENOME_DIR / ".bowtie1_index_complete"),
    params:
        prefix=config["reference"]["bowtie1_index"],
    log:
        str(LOG_DIR / "align" / "bowtie1_index.log"),
    shell:
        """
        bowtie-build --threads {threads} {input.fa} {params.prefix} > {log} 2>&1
        touch {output.done}
        """


rule three_letter_references:
    input:
        fa=config["reference"]["genome_fa"],
    output:
        c2t=str(THREE_LETTER_C2T_FA),
        g2a=str(THREE_LETTER_G2A_FA),
        done=str(THREE_LETTER_DIR / ".references_complete"),
    params:
        script=str(SRNATAPS_ROOT / "three_letter.py"),
        outdir=str(THREE_LETTER_DIR),
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} convert-reference \
            --input {input.fa} \
            --c2t {output.c2t}.tmp \
            --g2a {output.g2a}.tmp
        mv {output.c2t}.tmp {output.c2t}
        mv {output.g2a}.tmp {output.g2a}
        touch {output.done}
        """


rule three_letter_c2t_index:
    threads:
        config.get("alignment", {}).get("threads", 8),
    input:
        fa=str(THREE_LETTER_C2T_FA),
        references=str(THREE_LETTER_DIR / ".references_complete"),
    output:
        done=str(THREE_LETTER_DIR / ".c2t_index_complete"),
    params:
        prefix=str(THREE_LETTER_C2T_INDEX),
    log:
        str(LOG_DIR / "three_letter" / "build_C2T.log"),
    resources:
        mem_mb=40000,
        runtime=480,
    shell:
        """
        mkdir -p {LOG_DIR}/three_letter
        bowtie-build --threads {threads} {input.fa} {params.prefix} \
            > {log} 2>&1
        touch {output.done}
        """


rule three_letter_g2a_index:
    threads:
        config.get("alignment", {}).get("threads", 8),
    input:
        fa=str(THREE_LETTER_G2A_FA),
        references=str(THREE_LETTER_DIR / ".references_complete"),
    output:
        done=str(THREE_LETTER_DIR / ".g2a_index_complete"),
    params:
        prefix=str(THREE_LETTER_G2A_INDEX),
    log:
        str(LOG_DIR / "three_letter" / "build_G2A.log"),
    resources:
        mem_mb=40000,
        runtime=480,
    shell:
        """
        mkdir -p {LOG_DIR}/three_letter
        bowtie-build --threads {threads} {input.fa} {params.prefix} \
            > {log} 2>&1
        touch {output.done}
        """


rule bowtie1_align:
    """
    Align short TAPS RNA reads with Bowtie1 seed mode.

    Short seeds allow candidate discovery when genuine TAPS conversions occur
    elsewhere in the read. Both genomic strands are searched because annotated
    RNAs occur on both strands. Independent Bowtie processes are used because
    current Bioconda Bowtie1 builds may ignore the -p thread setting.
    """
    threads:
        config.get("alignment", {}).get("threads", 8),
    input:
        fq=str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        done=alignment_index_markers,
    output:
        bam=str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai=str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
    params:
        strategy=ALIGNMENT_STRATEGY,
        index=config["reference"]["bowtie1_index"],
        c2t_index=str(THREE_LETTER_C2T_INDEX),
        g2a_index=str(THREE_LETTER_G2A_INDEX),
        three_letter_script=str(SRNATAPS_ROOT / "three_letter.py"),
        three_letter_runner=str(SRNATAPS_ROOT / "workflow" / "scripts" / "three_letter_align.sh"),
        alignment_args=bowtie_alignment_args(),
        strand_arg=bowtie_strand_arg(),
        multimappers=config.get("alignment", {}).get("multimappers", 1),
        best_strata=(
            "--best --strata"
            if config.get("alignment", {}).get("best_strata", False)
            else ""
        ),
        max_multi_arg=(
            f"-m {config.get('alignment', {}).get('max_multimappers')}"
            if config.get("alignment", {}).get("max_multimappers")
            else ""
        ),
    resources:
        mem_mb=lambda wc, input: est_mem(
            8000, input.fq, scale=8, floor_mb=40000, ceil_mb=64000
        ),
        runtime=480,
    log:
        bowtie=str(LOG_DIR / "align" / "{sample}_bowtie.log"),
        flagstat=str(LOG_DIR / "align" / "{sample}_flagstat.log"),
        plus=str(LOG_DIR / "three_letter" / "{sample}_plus.log"),
        minus=str(LOG_DIR / "three_letter" / "{sample}_minus.log"),
    shell:
        """
        if [ "{params.strategy}" = "three_letter" ]; then
            bash {params.three_letter_runner} \
                {input.fq} \
                {output.bam} \
                {params.c2t_index} \
                {params.g2a_index} \
                "$(command -v python)" \
                {params.three_letter_script} \
                {threads} \
                {log.bowtie} \
                {log.flagstat} \
                {log.plus} \
                {log.minus} \
                "$TMPDIR"
            exit 0
        fi

        mkdir -p {ALIGN_DIR}
        WORK=$(mktemp -d "$TMPDIR/srnataps.{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$WORK"' EXIT

        gzip -dc {input.fq} \
        | awk -v out="$WORK" -v n={threads} '
            {{
                chunk = int((NR - 1) / 4) % n
                file = sprintf("%s/chunk_%03d.fastq", out, chunk)
                print >> file
            }}
        '

        PIDS=""
        for fq in "$WORK"/chunk_*.fastq; do
            (
                bowtie \
                    {params.alignment_args} \
                    {params.strand_arg} \
                    -k {params.multimappers} \
                    {params.best_strata} \
                    {params.max_multi_arg} \
                    -p 1 \
                    -q \
                    --sam \
                    -x {params.index} \
                    "$fq" \
                    2> "$fq.bowtie.log" \
                | samtools view -bS -F 4 -o "$fq.bam" -
            ) &
            PIDS="$PIDS $!"
        done

        STATUS=0
        for pid in $PIDS; do
            wait "$pid" || STATUS=1
        done
        if [ "$STATUS" -ne 0 ]; then
            cat "$WORK"/chunk_*.bowtie.log > {log.bowtie}
            echo "One or more Bowtie workers failed" >&2
            exit "$STATUS"
        fi

        awk '
            /# reads processed:/ {{
                value = $0; sub(/^.*: /, "", value); sub(/ .*/, "", value)
                total += value
            }}
            /# reads with at least one (reported )?alignment:/ {{
                value = $0; sub(/^.*: /, "", value); sub(/ .*/, "", value)
                mapped += value
            }}
            /^Reported / {{ reported += $2 }}
            END {{
                failed = total - mapped
                mapped_pct = total ? 100 * mapped / total : 0
                failed_pct = total ? 100 * failed / total : 0
                print "# reads processed: " total
                printf "# reads with at least one reported alignment: %d (%.2f%%)", mapped, mapped_pct
                print ""
                printf "# reads that failed to align: %d (%.2f%%)", failed, failed_pct
                print ""
                print "Reported " reported " alignments"
            }}
        ' "$WORK"/chunk_*.bowtie.log > {log.bowtie}

        samtools cat -o "$WORK/combined.bam" "$WORK"/chunk_*.fastq.bam
        samtools sort -@ 1 -o {output.bam} "$WORK/combined.bam"
        samtools index {output.bam} {output.bai}
        samtools flagstat {output.bam} > {log.flagstat}
        """
