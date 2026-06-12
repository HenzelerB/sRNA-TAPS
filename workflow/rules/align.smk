# =============================================================================
# rules/align.smk — Genome setup + Bowtie1 index + alignment
#
# For new users, the pipeline downloads and indexes the genome automatically.
# Existing users with a genome already in place: set the paths in config.yaml
# and the download steps will be skipped (Snakemake sees outputs already exist).
#
# Ensembl release is set in config["reference"]["ensembl_release"] (default 112).
# =============================================================================

ENSEMBL_RELEASE = config["reference"].get("ensembl_release", 112)
ENSEMBL_FA_URL  = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/"
    f"Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
)
ENSEMBL_GTF_URL = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf/homo_sapiens/"
    f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
)


rule download_genome_fa:
    """
    Download Ensembl GRCh38 unmasked FASTA.
    Unmasked (dna.toplevel) is required for TAPS — the repeat-masked version
    soft-masks cytosines, which would suppress genuine C→T modification signal.
    Skipped automatically if the FASTA already exists.
    """
    output:
        fa = config["reference"]["genome_fa"],
    params:
        url    = ENSEMBL_FA_URL,
        fa_gz  = config["reference"]["genome_fa"] + ".gz",
        outdir = str(GENOME_DIR),
    log:
        str(LOG_DIR / "align" / "download_genome_fa.log"),
    shell:
        """
        mkdir -p {params.outdir}
        echo "Downloading genome FASTA from Ensembl release {ENSEMBL_RELEASE}..." > {log}
        wget -q --show-progress -O {params.fa_gz} {params.url} >> {log} 2>&1
        echo "Decompressing..." >> {log}
        gunzip -f {params.fa_gz} >> {log} 2>&1
        echo "Done: {output.fa}" >> {log}
        """


rule download_gtf:
    """
    Download Ensembl GRCh38 GTF annotation.
    Used for biotype-split annotation.
    Skipped automatically if the GTF already exists.
    """
    output:
        gtf = config["reference"]["gtf"],
    params:
        url     = ENSEMBL_GTF_URL,
        gtf_gz  = config["reference"]["gtf"] + ".gz",
        outdir  = str(GENOME_DIR),
    log:
        str(LOG_DIR / "align" / "download_gtf.log"),
    shell:
        """
        mkdir -p {params.outdir}
        echo "Downloading GTF from Ensembl release {ENSEMBL_RELEASE}..." > {log}
        wget -q --show-progress -O {params.gtf_gz} {params.url} >> {log} 2>&1
        echo "Decompressing..." >> {log}
        gunzip -f {params.gtf_gz} >> {log} 2>&1
        echo "Done: {output.gtf}" >> {log}
        """


rule index_genome_fa:
    """
    samtools faidx — required by caller.py for random FASTA access.
    """
    input:
        fa = config["reference"]["genome_fa"],
    output:
        fai = config["reference"]["genome_fa"] + ".fai",
    log:
        str(LOG_DIR / "align" / "index_genome_fa.log"),
    shell:
        """
        samtools faidx {input.fa} > {log} 2>&1
        """


rule bowtie1_index:
    """
    Build Bowtie1 genome index — run once.
    Unmasked FASTA is required (see download_genome_fa rationale).
    """
    input:
        fa = config["reference"]["genome_fa"],
    output:
        done = str(GENOME_DIR / ".bowtie1_index_complete"),
    params:
        prefix = config["reference"]["bowtie1_index"],
    log:
        str(LOG_DIR / "align" / "bowtie1_index.log"),
    shell:
        """
        bowtie-build --threads {threads} {input.fa} {params.prefix} > {log} 2>&1
        touch {output.done}
        """


rule bowtie1_align:
    """
    Bowtie1 alignment — TAPS-aware, small RNA.
    -v 2            : allow 2 mismatches (tolerates genuine C→T TAPS signal)
    --norc          : strand-specific TruSeq small RNA libraries
    -k 10           : report up to 10 alignments (tRNA multi-mappers + XA tag)
    --best --strata : report only best-stratum hits
    -m 100          : discard reads with >100 valid alignments (repetitive elements)
    """
    input:
        fq   = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        done = str(GENOME_DIR / ".bowtie1_index_complete"),
    output:
        bam = str(ALIGN_DIR / "{sample}.sorted.bam"),
        bai = str(ALIGN_DIR / "{sample}.sorted.bam.bai"),
    params:
        index = config["reference"]["bowtie1_index"],
    resources:
        mem_mb  = lambda wc, input: est_mem(
            8000, input.fq,
            scale=8, floor_mb=12000, ceil_mb=64000
        ),
        runtime = 2880,
    log:
        bowtie   = str(LOG_DIR / "align" / "{sample}_bowtie.log"),
        flagstat = str(LOG_DIR / "align" / "{sample}_flagstat.log"),
    shell:
        """
        mkdir -p {ALIGN_DIR}
        bowtie \
            -v 2 \
            --norc \
            -k 10 \
            --best \
            --strata \
            -m 100 \
            -p {threads} \
            --sam \
            {params.index} \
            {input.fq} \
            2> {log.bowtie} \
        | samtools view -bS -F 4 - \
        | samtools sort -@ 4 -o {output.bam}

        samtools index {output.bam} {output.bai}
        samtools flagstat {output.bam} > {log.flagstat}
        """
