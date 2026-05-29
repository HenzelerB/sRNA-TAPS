# =============================================================================
# rules/align.smk — Bowtie1 genome index and alignment
# =============================================================================

rule bowtie1_index:
    """
    Build Bowtie1 genome index from reference FASTA.

    Why Bowtie1 for small RNA:
        - Reads are 18-50 nt — too short for bwa-mem2 minimum seed length
        - -v mode tolerates total mismatches without quality weighting,
          correctly handling TAPS C→T at modified sites as real signal
        - Splice-awareness (STAR/HISAT2) adds no value for non-spliced small RNA
    """
    input:
        fa = config["reference"]["genome_fa"],
    output:
        done  = str(GENOME_DIR / ".bowtie1_index_complete"),
        index = str(GENOME_DIR / "genome.1.ebwt"),
    params:
        outdir = str(GENOME_DIR),
        prefix = str(GENOME_DIR / "genome"),
    log:
        str(LOG_DIR / "align" / "bowtie1_index.log"),
    threads: 20
    resources:
        mem_mb   = 64000,
        runtime  = 720,
    shell:
        """
        mkdir -p {params.outdir}
        bowtie-build --threads {threads} {input.fa} {params.prefix} > {log} 2>&1
        touch {output.done}
        """


rule bowtie1_align:
    """
    Bowtie1 alignment of trimmed small RNA reads.

    Parameter rationale:
        -v 2         : up to 2 total mismatches (TAPS C→T are genuine mismatches)
        --norc       : strand-specific (TruSeq small RNA = sense strand only)
                       prevents false G→A calls on reverse complement
        -k 10        : report up to 10 alignments (tRNA ~600 near-identical copies)
        --best       : report alignments in the best-scoring stratum only
        --strata     : prevents mixing 0- and 2-mismatch hits in XA count
        -m 100       : discard reads mapping to >100 loci (uninterpretable repeats)
        --sam        : SAM output required for XA:i:N tag (multi-mapper weighting)
    """
    input:
        fq    = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        index = str(GENOME_DIR / ".bowtie1_index_complete"),
    output:
        bam   = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai   = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
    params:
        prefix    = str(GENOME_DIR / "genome"),
        mm        = config["alignment"]["mismatches"],
        multi     = config["alignment"]["multimappers"],
        max_multi = config["alignment"]["max_multimappers"],
    log:
        bowtie  = str(LOG_DIR / "align" / "{sample}_bowtie.log"),
        flagstat= str(LOG_DIR / "align" / "{sample}_flagstat.log"),
    threads: 8
    resources:
        mem_mb   = 40000,
        runtime  = 480,
    shell:
        """
        mkdir -p {ALIGN_DIR}
        bowtie \
            -v {params.mm} \
            --norc \
            -k {params.multi} \
            --best \
            --strata \
            -m {params.max_multi} \
            -p {threads} \
            --sam \
            {params.prefix} \
            {input.fq} \
            2> {log.bowtie} \
        | samtools view -bS -F 4 - \
        | samtools sort -@ 4 -o {output.bam}

        samtools index {output.bam}
        samtools flagstat {output.bam} > {log.flagstat}
        """
