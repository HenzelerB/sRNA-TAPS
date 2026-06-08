# =============================================================================
# rules/compare.smk — Benchmarking comparison
# =============================================================================

rule compare:
    """
    Compare custom TAPS pipeline output vs rastair, asTair, Bismark.
    Produces concordance and correlation summary tables.
    Bismark rates are inverted (1 - rate) before comparison.
    """
    input:
        calls   = expand(str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
                         sample=SAMPLE_NAMES, biotype=BIOTYPES),
        rastair = expand(str(BENCH_DIR / "rastair" / "all" / "{sample}" / "{sample}.bed.gz"),
                         sample=SAMPLE_NAMES),
        astair  = expand(str(BENCH_DIR / "astair" / "{biotype}" / ".done_{sample}"),
                         sample=SAMPLE_NAMES, biotype=BIOTYPES),
        bismark = expand(str(BENCH_DIR / "bismark" / "{sample}" / ".done"),
                         sample=SAMPLE_NAMES),
    output:
        concordance  = str(COMPARE_DIR / "concordance_summary.tsv"),
        correlation  = str(COMPARE_DIR / "correlation_summary.tsv"),
    params:
        script      = str(Path(workflow.basedir).parent / "srnataps" / "compare.py"),
        calls_dir   = str(CALLS_DIR),
        bench_dir   = str(BENCH_DIR),
        out_dir     = str(COMPARE_DIR),
    log:
        str(LOG_DIR / "compare" / "compare.log"),
    shell:
        """
        mkdir -p {params.out_dir}
        python {params.script} \
            --calls-dir  {params.calls_dir} \
            --bench-dir  {params.bench_dir} \
            --out-dir    {params.out_dir} \
            > {log} 2>&1
        """
