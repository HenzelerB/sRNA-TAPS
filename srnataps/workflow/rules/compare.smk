# =============================================================================
# rules/compare.smk — Benchmarking comparison: sRNA-TAPS vs rastair/asTair/Bismark
# =============================================================================

rule compare:
    """
    Compare sRNA-TAPS output vs rastair, asTair, Bismark.
    Produces concordance and correlation summary tables.
    Bismark rates are inverted (1 - rate) before comparison.
    """
    input:
        calls   = expand(str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
                         sample=SAMPLES, biotype=BIOTYPES),
        rastair = expand(str(BENCH_DIR / "rastair" / "all" / "{sample}" / "{sample}.bed.gz"),
                         sample=SAMPLES),
        rastair_cpg = expand(str(BENCH_DIR / "rastair" / "cpg" / "{sample}" / "{sample}.bed.gz"),
                         sample=SAMPLES),
        astair  = expand(str(BENCH_DIR / "astair" / "{biotype}" / ".done_{sample}"),
                         sample=SAMPLES, biotype=BIOTYPES),
        bismark = expand(str(BENCH_DIR / "bismark" / "{sample}" / ".done"),
                         sample=SAMPLES),
    output:
        concordance = str(COMPARE_DIR / "concordance_summary.tsv"),
        correlation = str(COMPARE_DIR / "correlation_summary.tsv"),
    params:
        script    = str(SRNATAPS_SCRIPTS / "compare.py"),
        calls_dir = str(CALLS_DIR),
        bench_dir = str(BENCH_DIR),
        out_dir   = str(COMPARE_DIR),
    log:
        str(LOG_DIR / "compare" / "compare.log"),
    shell:
        """
        mkdir -p {params.out_dir}
        cd $(dirname {params.out_dir})
        python {params.script} > {log} 2>&1
        """
