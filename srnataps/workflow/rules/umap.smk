# =============================================================================
# rules/umap.smk — UMAP dimensionality reduction and cluster analysis
# =============================================================================
# Generates comprehensive visualizations of sample relationships based on
# modification patterns across biotypes and conditions.
#
# Position: 07a (immediate QC after base taps_calls)
#
# Outputs:
# - umap_coordinates.tsv: 2D UMAP embedding with cluster assignments
# - outlier_classification.tsv: Outlier detection results
# - neighbors.json: Sample similarity network (k-NN)
# - feature_matrix.tsv: Engineered features used for analysis
# - *.png: Visualizations (clusters, outliers, neighborhoods, metrics)
# =============================================================================

# Define UMAP script location
SRNATAPS_SCRIPTS_UMAP = Path(workflow.basedir).parent

rule umap_analysis:
    """
    Generate UMAP analysis of modification patterns across all samples.

    Position: Early QC stage (07a) after individual sample calls.
    Uses base per-sample modification calls to create feature matrix,
    then applies UMAP for dimensionality reduction and comprehensive
    cluster/outlier/neighbor analysis.

    This is a quality control gate: verifies that samples separate by
    condition/cell-line using raw modification signal before condition analysis.
    """
    input:
        expand(
            str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
            sample=SAMPLES, biotype=BIOTYPES,
        ),
    output:
        coords = str(UMAP_DIR / "umap_coordinates.tsv"),
        dbscan_membership = str(UMAP_DIR / "dbscan_membership.tsv"),
        outliers = str(UMAP_DIR / "outlier_classification.tsv"),
        neighbors = str(UMAP_DIR / "neighbors.json"),
        features = str(UMAP_DIR / "feature_matrix.tsv"),
        condition_statistics = str(UMAP_DIR / "condition_statistics.tsv"),
        pairwise_statistics = str(UMAP_DIR / "pairwise_statistics.tsv"),
        global_statistics = str(UMAP_DIR / "global_statistics.json"),
        plot_clusters = str(UMAP_DIR / "01_umap_clusters.png"),
        plot_clusters_pdf = str(UMAP_DIR / "01_umap_clusters.pdf"),
        plot_outliers = str(UMAP_DIR / "02_umap_outliers.png"),
        plot_neighbors = str(UMAP_DIR / "03_umap_neighbors.png"),
        plot_metrics = str(UMAP_DIR / "04_umap_metrics.png"),
    params:
        script = str(SRNATAPS_SCRIPTS_UMAP / "umap_analysis.py"),
        calls_dir = str(CALLS_DIR),
        biotypes = ",".join(BIOTYPES),
        samples_tsv = "samples.tsv" if Path("samples.tsv").exists() else "",
        outdir = str(UMAP_DIR),
        n_neighbors = config.get("umap", {}).get("n_neighbors", 15),
        eps = config.get("umap", {}).get("eps", 0.5),
        random_state = config.get("umap", {}).get("random_state", 42),
        method = config.get("umap", {}).get("method", "umap"),
    resources:
        mem_mb = 8000,
        runtime = 240,
    log:
        str(LOG_DIR / "umap" / "umap_analysis.log"),
    shell:
        """
        mkdir -p $(dirname {output.coords}) $(dirname {log})

        # Check if we have the required packages
        if ! python -c "import umap, pacmap, matplotlib, seaborn, sklearn, scipy" 2>/dev/null; then
            echo "Installing required packages for UMAP analysis..."
            pip install umap-learn pacmap matplotlib seaborn scikit-learn scipy 2>&1 | tail -5
        fi

        python {params.script} \
            --calls-dir {params.calls_dir} \
            --biotypes {params.biotypes} \
            --outdir {params.outdir} \
            --n-neighbors {params.n_neighbors} \
            --eps {params.eps} \
            --random-state {params.random_state} \
            --method {params.method} \
            $([ -n "{params.samples_tsv}" ] && echo "--samples-tsv {params.samples_tsv}") \
            > {log} 2>&1

        if [ $? -ne 0 ]; then
            echo "UMAP analysis failed. Check log for details: {log}"
            exit 1
        fi
        """


# Optional: UMAP on stringent condition calls for validation
rule umap_analysis_stringent:
    """
    Alternative UMAP using stringent condition-level calls (07e).
    Useful for validating that statistical filtering preserves sample structure.
    Produces outputs with _stringent suffix.
    """
    input:
        expand(
            str(STRINGENT_DIR / "{biotype}" / "treat_{cell_line}_stringent_{biotype}_taps.tsv"),
            cell_line=CELL_LINES, biotype=BIOTYPES,
        ) if CONDITION_ANALYSIS_ENABLED else [],
    output:
        coords = str(UMAP_DIR / "umap_coordinates_stringent.tsv"),
        plot_clusters = str(UMAP_DIR / "umap_clusters_stringent.png"),
    params:
        script = str(SRNATAPS_SCRIPTS_UMAP / "umap_analysis.py"),
        stringent_dir = str(STRINGENT_DIR),
        biotypes = ",".join(BIOTYPES),
        outdir_tmp = "/tmp/umap_stringent",
        outdir = str(UMAP_DIR),
    resources:
        mem_mb = 8000,
        runtime = 180,
    log:
        str(LOG_DIR / "umap" / "umap_stringent.log"),
    shell:
        """
        mkdir -p {params.outdir_tmp}

        python {params.script} \
            --stringent-dir {params.stringent_dir} \
            --biotypes {params.biotypes} \
            --outdir {params.outdir_tmp} \
            > {log} 2>&1

        mv {params.outdir_tmp}/umap_coordinates.tsv {output.coords}
        mv {params.outdir_tmp}/01_umap_clusters.png {params.outdir}/umap_clusters_stringent.png
        """


rule umap_analysis_pb_tet:
    """Compare PB+ TET- and PB+ TET+ samples without untreated samples."""
    input:
        expand(
            str(CALLS_DIR / "{biotype}" / "{sample}_{biotype}_taps.tsv"),
            sample=SAMPLES, biotype=BIOTYPES,
        ),
    output:
        coords=str(UMAP_DIR / "pb_tet_comparison" / "umap_coordinates.tsv"),
        dbscan_membership=str(UMAP_DIR / "pb_tet_comparison" / "dbscan_membership.tsv"),
        outliers=str(UMAP_DIR / "pb_tet_comparison" / "outlier_classification.tsv"),
        neighbors=str(UMAP_DIR / "pb_tet_comparison" / "neighbors.json"),
        features=str(UMAP_DIR / "pb_tet_comparison" / "feature_matrix.tsv"),
        condition_statistics=str(UMAP_DIR / "pb_tet_comparison" / "condition_statistics.tsv"),
        pairwise_statistics=str(UMAP_DIR / "pb_tet_comparison" / "pairwise_statistics.tsv"),
        global_statistics=str(UMAP_DIR / "pb_tet_comparison" / "global_statistics.json"),
        plot_clusters=str(UMAP_DIR / "pb_tet_comparison" / "01_umap_clusters.png"),
        plot_clusters_pdf=str(UMAP_DIR / "pb_tet_comparison" / "01_umap_clusters.pdf"),
        plot_outliers=str(UMAP_DIR / "pb_tet_comparison" / "02_umap_outliers.png"),
        plot_neighbors=str(UMAP_DIR / "pb_tet_comparison" / "03_umap_neighbors.png"),
        plot_metrics=str(UMAP_DIR / "pb_tet_comparison" / "04_umap_metrics.png"),
    params:
        script=str(SRNATAPS_SCRIPTS_UMAP / "umap_analysis.py"),
        calls_dir=str(CALLS_DIR),
        biotypes=",".join(BIOTYPES),
        samples_tsv="samples.tsv" if Path("samples.tsv").exists() else "",
        outdir=str(UMAP_DIR / "pb_tet_comparison"),
        n_neighbors=config.get("umap", {}).get("n_neighbors", 15),
        eps=config.get("umap", {}).get("eps", 0.5),
        random_state=config.get("umap", {}).get("random_state", 42),
        method=config.get("umap", {}).get("method", "umap"),
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        str(LOG_DIR / "umap" / "umap_pb_tet.log"),
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir} $(dirname {log})
        python {params.script} \
            --calls-dir {params.calls_dir} \
            --biotypes {params.biotypes} \
            --outdir {params.outdir} \
            --n-neighbors {params.n_neighbors} \
            --eps {params.eps} \
            --random-state {params.random_state} \
            --method {params.method} \
            --conditions pb_ctrl,treat \
            $([ -n "{params.samples_tsv}" ] && echo "--samples-tsv {params.samples_tsv}") \
            > {log} 2>&1
        """
