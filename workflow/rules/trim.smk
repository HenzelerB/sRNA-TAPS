# =============================================================================
# rules/trim.smk — Trim Galore adapter trimming
# =============================================================================

rule trim_galore:
    input:
        fq = lambda wc: SAMPLES_DF.loc[wc.sample, "fastq"],
    output:
        fq  = str(TRIM_DIR / "{sample}_trimmed.fq.gz"),
        log = str(TRIM_DIR / "{sample}_trim_galore.log"),
    log:
        str(LOG_DIR / "trim" / "{sample}_trim_galore.log"),
    params:
        adapter    = config["trimming"]["adapter"],
        min_length = config["trimming"]["min_length"],
    shell:
        """
        mkdir -p {TRIM_DIR}
        # Default TruSeq small-RNA adapter -> use --small_rna (preserves the
        # small-RNA-specific short-read removal behaviour). A custom adapter
        # set in config.yaml -> pass it explicitly with -a instead.
        DEFAULT_SMALLRNA_ADAPTER="TGGAATTCTCGGGTGCCAAGG"
        if [ "{params.adapter}" = "$DEFAULT_SMALLRNA_ADAPTER" ]; then
            ADAPTER_ARG="--small_rna"
        else
            ADAPTER_ARG="-a {params.adapter}"
        fi
        trim_galore \
            $ADAPTER_ARG \
            --length {params.min_length} \
            --cores {threads} \
            -o {TRIM_DIR} \
            {input.fq} \
            > {log} 2>&1

        # Normalise output name: trim_galore adds _trimmed suffix to the stem
        STEM=$(basename {input.fq} .fq.gz | sed 's/\.merged$//')
        SRC="{TRIM_DIR}/${{STEM}}_trimmed.fq.gz"
        DST="{output.fq}"
        [ -f "$SRC" ] && [ "$SRC" != "$DST" ] && mv "$SRC" "$DST" || true

        cp {log} {output.log}
        """
