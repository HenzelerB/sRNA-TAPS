#!/usr/bin/env python3
"""
UMAP analysis for sRNA-TAPS modifications.

Generates comprehensive visualizations showing:
- Clusters: sample grouping by modification patterns
- Outliers: samples with unusual modification profiles
- Local neighbors: sample similarity networks
- Global layout: high-dimensional structure reduction
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
import json
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from itertools import combinations

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

AURORA_CONDITION_COLORS = {
    "no_treat": "#4DDEB8",
    "pb_ctrl": "#5BAFD0",
    "treat": "#FFD680",
}
CONDITION_LABELS = {
    "no_treat": "PB- TET-",
    "pb_ctrl": "PB+ TET-",
    "treat": "PB+ TET+",
}
CONDITION_ORDER = ("no_treat", "pb_ctrl", "treat")
PLOT_ALPHA = 0.6

try:
    import umap
    import pacmap
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import zscore, kruskal, mannwhitneyu
except ImportError as e:
    logger.error(f"Missing required package: {e}")
    logger.error("Install: pip install umap-learn matplotlib seaborn scikit-learn scipy")
    sys.exit(1)


def load_calls(calls_dir: Path, biotypes: List[str], cell_lines: List[str]) -> pd.DataFrame:
    """
    Load modification calls across all samples and biotypes.
    Works with both individual sample calls (07.taps_calls) and condition-level calls (07e.stringent).

    Returns DataFrame with rows=samples, columns=features (modification rates by biotype).
    """
    logger.info(f"Loading calls from {calls_dir}")

    # Collect modification rates per sample per biotype
    data_dict = {}

    for biotype in biotypes:
        biotype_dir = calls_dir / biotype
        if not biotype_dir.exists():
            logger.warning(f"Biotype directory not found: {biotype_dir}")
            continue

        # Find all call files
        call_files = list(biotype_dir.glob("*_taps.tsv"))

        for call_file in call_files:
            # Parse filename: various formats depending on stage
            # Individual: {sample}_{biotype}_taps.tsv
            # Condition: treat_{CELLLINE}_stringent_{biotype}_taps.tsv
            parts = call_file.stem.split("_")

            try:
                df = pd.read_csv(call_file, sep="\t", low_memory=False)
                if df.empty:
                    logger.debug(f"Empty file: {call_file.name}")
                    continue

                # Calculate modification rate (fraction of called sites with modifications)
                if "mod_rate" in df.columns:
                    mod_rate = df["mod_rate"].mean()
                else:
                    # Alternative: use binary calls
                    mod_rate = (df.iloc[:, -1] > 0).sum() / len(df) if len(df) > 0 else 0

                # Extract sample/condition name from filename
                if call_file.name.endswith("_taps.tsv"):
                    sample_key = call_file.stem.replace(f"_{biotype}_taps", "")
                else:
                    sample_key = "_".join(parts[:-2]) if len(parts) > 2 else parts[0]

                if sample_key not in data_dict:
                    data_dict[sample_key] = {}

                feature_key = f"{biotype}_mod_rate"
                data_dict[sample_key][feature_key] = mod_rate

                logger.debug(f"Loaded {call_file.name}: {biotype} mod_rate={mod_rate:.4f}")

            except Exception as e:
                logger.warning(f"Failed to load {call_file.name}: {e}")
                continue

    # Convert to DataFrame
    if not data_dict:
        logger.error("No calls loaded!")
        return pd.DataFrame()

    feature_df = pd.DataFrame.from_dict(data_dict, orient="index")
    logger.info(f"Loaded features for {len(feature_df)} samples x {len(feature_df.columns)} biotypes")

    return feature_df


def engineer_features(feature_df: pd.DataFrame,
                     calls_dir: Path,
                     biotypes: List[str],
                     samples_tsv: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Engineer additional features from raw calls:
    - Coverage per biotype
    - Modification density
    - Cross-biotype correlation patterns
    - Sample metadata
    """
    logger.info("Engineering additional features...")

    if feature_df.empty:
        return feature_df, pd.DataFrame()

    extended_features = feature_df.copy()

    # Load sample metadata if available
    metadata = None
    if samples_tsv and samples_tsv.exists():
        metadata = pd.read_csv(samples_tsv, sep="\t")
        metadata = metadata.set_index("sample")

    # For each biotype, add coverage-weighted features
    for biotype in biotypes:
        biotype_dir = calls_dir / biotype
        if not biotype_dir.exists():
            continue

        coverage_data = {}
        for call_file in biotype_dir.glob("*_taps.tsv"):
            try:
                df = pd.read_csv(call_file, sep="\t", low_memory=False)
                if not df.empty and "coverage" in df.columns:
                    sample_key = "_".join(call_file.stem.split("_")[:-3])
                    coverage_data[sample_key] = df["coverage"].mean()
            except:
                continue

        if coverage_data:
            cov_series = pd.Series(coverage_data)
            extended_features[f"{biotype}_avg_coverage"] = cov_series

    # Fill missing values with 0
    extended_features = extended_features.fillna(0)

    return extended_features, metadata


def select_replicates_per_group(features: pd.DataFrame, metadata: pd.DataFrame,
                                replicates_per_group: int) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Retain the profiles closest to each condition-by-cell-line centroid."""
    required = {"condition", "cell_line"}
    if metadata is None or not required.issubset(metadata.columns):
        raise ValueError("Replicate selection requires condition and cell_line metadata")
    metadata = metadata.reindex(features.index)
    if metadata[list(required)].isna().any().any():
        raise ValueError("Replicate selection found samples without condition or cell-line metadata")
    scaled = pd.DataFrame(StandardScaler().fit_transform(features),
                          index=features.index, columns=features.columns)
    selection_rows = []
    retained = []
    for (condition, cell_line), group_metadata in metadata.groupby(["condition", "cell_line"], sort=False):
        group_samples = group_metadata.index.tolist()
        if len(group_samples) < replicates_per_group:
            raise ValueError(
                f"{condition}/{cell_line} has {len(group_samples)} replicates; "
                f"cannot retain {replicates_per_group}"
            )
        centroid = scaled.loc[group_samples].mean(axis=0)
        distances = ((scaled.loc[group_samples] - centroid) ** 2).sum(axis=1) ** 0.5
        retained_group = distances.sort_values(kind="stable").head(replicates_per_group).index
        retained.extend(retained_group)
        for sample, distance in distances.items():
            selection_rows.append({
                "sample": sample,
                "condition": condition,
                "condition_label": CONDITION_LABELS.get(condition, condition),
                "cell_line": cell_line,
                "centroid_distance": distance,
                "selection": "retained" if sample in retained_group else "excluded",
            })
    selection_df = pd.DataFrame(selection_rows).sort_values(
        ["condition", "cell_line", "selection", "centroid_distance"]
    )
    return features.loc[retained], metadata.loc[retained], selection_df


def compute_umap(features: pd.DataFrame, n_neighbors: int = 15,
                min_dist: float = 0.1, n_components: int = 2,
                random_state: int = 42, method: str = "umap") -> np.ndarray:
    """
    Compute UMAP embedding with specified parameters.

    Returns: (n_samples, n_components) array of UMAP coordinates
    """
    logger.info(f"Computing embedding with method={method}, n_neighbors={n_neighbors}")

    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)

    if method == "pca":
        reducer = PCA(n_components=n_components, random_state=random_state)
        embedding = reducer.fit_transform(features_scaled)
        logger.info("PCA explained variance: %s", reducer.explained_variance_ratio_)
        return embedding
    if method == "pacmap":
        reducer = pacmap.PaCMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            MN_ratio=0.5,
            FP_ratio=2.0,
            random_state=random_state,
        )
    elif method == "umap":
        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric="euclidean",
            random_state=random_state,
            verbose=1,
        )
    else:
        raise ValueError("method must be 'umap' or 'pacmap'")

    embedding = reducer.fit_transform(features_scaled)
    logger.info(f"UMAP computed: shape={embedding.shape}")

    return embedding


def detect_clusters(embedding: np.ndarray, eps: float = 0.5,
                   min_samples: int = 3) -> np.ndarray:
    """
    Detect clusters using DBSCAN on UMAP embedding.

    Returns: cluster labels (-1 = noise/outlier)
    """
    logger.info(f"Detecting clusters with DBSCAN (eps={eps}, min_samples={min_samples})")

    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    labels = dbscan.fit_predict(embedding)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_outliers = list(labels).count(-1)

    logger.info(f"Found {n_clusters} clusters with {n_outliers} outliers")

    return labels


def identify_outliers(embedding: np.ndarray, features: pd.DataFrame,
                     zscore_threshold: float = 3.0) -> Dict[str, List]:
    """
    Identify outliers using multiple criteria:
    - Statistical outliers (high z-score in feature space)
    - Isolation in UMAP space (low local density)
    """
    outliers = {"statistical": [], "isolation": [], "combined": []}

    # Statistical outliers (z-score > threshold in any feature)
    feature_zscores = np.abs(zscore(features.fillna(0), axis=0))
    statistical_outlier_idx = np.where(feature_zscores.max(axis=1) > zscore_threshold)[0]
    outliers["statistical"] = features.index[statistical_outlier_idx].tolist()

    # Isolation outliers (sparse in UMAP space)
    neigh = NearestNeighbors(n_neighbors=5)
    neigh.fit(embedding)
    distances, _ = neigh.kneighbors(embedding)
    mean_dist = distances[:, -1]  # Distance to 5th nearest neighbor
    isolation_threshold = mean_dist.mean() + 2 * mean_dist.std()
    isolation_outlier_idx = np.where(mean_dist > isolation_threshold)[0]
    outliers["isolation"] = features.index[isolation_outlier_idx].tolist()

    # Combined outliers
    outliers["combined"] = list(set(outliers["statistical"] + outliers["isolation"]))

    logger.info(f"Identified outliers: statistical={len(outliers['statistical'])}, "
               f"isolation={len(outliers['isolation'])}, combined={len(outliers['combined'])}")

    return outliers


def compute_local_neighbors(embedding: np.ndarray, features: pd.DataFrame,
                           n_neighbors: int = 5) -> Dict[str, List]:
    """
    Compute k-nearest neighbors for each sample in UMAP space.

    Returns: dict mapping sample name -> list of neighbor sample names
    """
    logger.info(f"Computing {n_neighbors}-nearest neighbors...")

    neigh = NearestNeighbors(n_neighbors=n_neighbors+1)  # +1 to include self
    neigh.fit(embedding)
    _, indices = neigh.kneighbors(embedding)

    neighbors = {}
    for i, sample in enumerate(features.index):
        # Skip the first neighbor (self) and get the rest
        neighbor_samples = [features.index[idx] for idx in indices[i][1:]]
        neighbors[sample] = neighbor_samples

    return neighbors


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Return BH-adjusted p-values while preserving the input index."""
    values = pvalues.to_numpy(dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for rank in range(len(values) - 1, -1, -1):
        position = order[rank]
        running = min(running, values[position] * len(values) / (rank + 1))
        adjusted[position] = running
    return pd.Series(np.clip(adjusted, 0, 1), index=pvalues.index)


def run_statistics(features: pd.DataFrame, metadata: Optional[pd.DataFrame],
                   outdir: Path, random_state: int = 42,
                   permutations: int = 999) -> None:
    """Run condition-level tests on modification-rate features."""
    if metadata is None or "condition" not in metadata:
        logger.warning("Skipping statistics: condition metadata unavailable")
        return
    conditions = metadata.reindex(features.index)["condition"]
    valid = conditions.notna() & conditions.isin(AURORA_CONDITION_COLORS)
    features = features.loc[valid]
    conditions = conditions.loc[valid]
    groups = {condition: features.loc[conditions == condition]
              for condition in CONDITION_ORDER
              if (conditions == condition).sum() > 0}
    omnibus_rows = []
    pairwise_rows = []
    for feature in features.columns:
        samples = [group[feature].dropna().to_numpy() for group in groups.values()]
        if len(samples) >= 2:
            statistic, pvalue = kruskal(*samples)
            omnibus_rows.append({"feature": feature, "test": "kruskal_wallis",
                                 "statistic": statistic, "p_value": pvalue,
                                 "n": len(features)})
        for first, second in combinations(groups, 2):
            first_values = groups[first][feature].dropna()
            second_values = groups[second][feature].dropna()
            statistic, pvalue = mannwhitneyu(first_values, second_values,
                                             alternative="two-sided")
            pairwise_rows.append({"feature": feature, "condition_a": first,
                                  "condition_a_label": CONDITION_LABELS[first],
                                  "condition_b": second, "test": "mann_whitney_u",
                                  "condition_b_label": CONDITION_LABELS[second],
                                  "statistic": statistic, "p_value": pvalue,
                                  "n_a": len(first_values), "n_b": len(second_values)})
    omnibus = pd.DataFrame(omnibus_rows)
    pairwise = pd.DataFrame(pairwise_rows)
    if not omnibus.empty:
        omnibus["q_value"] = benjamini_hochberg(omnibus["p_value"])
    if not pairwise.empty:
        pairwise["q_value"] = benjamini_hochberg(pairwise["p_value"])
    omnibus.to_csv(outdir / "condition_statistics.tsv", sep="\t", index=False)
    pairwise.to_csv(outdir / "pairwise_statistics.tsv", sep="\t", index=False)

    rng = np.random.default_rng(random_state)
    matrix = StandardScaler().fit_transform(features.to_numpy(dtype=float))
    centroid_distance = 0.0
    overall = matrix.mean(axis=0)
    for group in groups.values():
        centroid_distance += len(group) * np.sum((group.to_numpy().mean(axis=0) - overall) ** 2)
    labels = conditions.to_numpy()
    permuted = np.zeros(permutations, dtype=float)
    for iteration in range(permutations):
        shuffled = rng.permutation(labels)
        value = 0.0
        for condition in groups:
            subset = matrix[shuffled == condition]
            if len(subset):
                value += len(subset) * np.sum((subset.mean(axis=0) - overall) ** 2)
        permuted[iteration] = value
    global_p = (1 + np.sum(permuted >= centroid_distance)) / (permutations + 1)
    with open(outdir / "global_statistics.json", "w") as handle:
        json.dump({"test": "permutation_between_condition_dispersion",
                   "statistic": float(centroid_distance),
                   "p_value": float(global_p), "permutations": permutations,
                   "n_samples": int(len(features)),
                   "conditions": {key: int(len(value)) for key, value in groups.items()}},
                  handle, indent=2)


def plot_umap_clusters(embedding: np.ndarray, features: pd.DataFrame,
                      cluster_labels: np.ndarray, outfile: Path,
                      metadata: Optional[pd.DataFrame] = None,
                      axis_labels: Tuple[str, str] = ("UMAP 1", "UMAP 2"),
                      title: str = "Sample-level modification profiles"):
    """Plot a publication-readable UMAP with condition and cell-line encoding."""
    logger.info(f"Plotting clusters -> {outfile}")
    metadata = metadata.reindex(features.index) if metadata is not None else None
    if metadata is not None and metadata.index.isna().any():
        logger.warning("Metadata missing for %d loaded samples", metadata.index.isna().sum())
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    cluster_ids = sorted(set(cluster_labels))
    cluster_palette = sns.color_palette("Greys", max(len(cluster_ids), 1) + 2)[2:]
    cluster_colors = {cluster: cluster_palette[i] for i, cluster in enumerate(cluster_ids)}
    conditions = ([condition for condition in CONDITION_ORDER
                   if condition in metadata["condition"].astype(str).unique()]
                  if metadata is not None and "condition" in metadata else ["all"])
    condition_palette = {
        condition: AURORA_CONDITION_COLORS.get(condition, "#1C4062")
        for condition in conditions
    }
    cell_lines = (metadata["cell_line"].astype(str).unique().tolist()
                  if metadata is not None and "cell_line" in metadata else ["all"])
    markers = dict(zip(cell_lines, ["o", "s", "^", "D", "P", "X"]))

    for i, sample in enumerate(features.index):
        condition = str(metadata.loc[sample, "condition"]) if metadata is not None and "condition" in metadata else "all"
        cell_line = str(metadata.loc[sample, "cell_line"]) if metadata is not None and "cell_line" in metadata else "all"
        cluster = int(cluster_labels[i])
        ax.scatter(embedding[i, 0], embedding[i, 1], s=72,
                   color=condition_palette[condition],
                   marker="x" if cluster == -1 else markers[cell_line],
                   linewidth=1.2 if cluster == -1 else 0,
                   alpha=PLOT_ALPHA,
                   zorder=3)

    condition_handles = [Line2D([0], [0], marker="o", color="w",
                                label=CONDITION_LABELS.get(condition, condition),
                                markerfacecolor=color, markersize=7)
                         for condition, color in condition_palette.items()]
    cell_handles = [Line2D([0], [0], marker=marker, color="0.25", label=cell_line,
                           linestyle="None", markersize=7)
                    for cell_line, marker in markers.items()]
    ax.legend(handles=condition_handles + cell_handles,
              title="Condition / cell line",
              bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.set_title(title)
    ax.grid(True, color="0.92", linewidth=0.7)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.savefig(outfile.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def plot_outliers(embedding: np.ndarray, features: pd.DataFrame,
                 outliers: Dict[str, List], outfile: Path,
                 axis_labels: Tuple[str, str] = ("UMAP 1", "UMAP 2")):
    """Plot UMAP highlighting outliers."""
    logger.info(f"Plotting outliers -> {outfile}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, (outlier_type, outlier_list) in enumerate(outliers.items()):
        if idx >= 2:
            break

        ax = axes[idx]

        # Color by outlier status
        colors = [AURORA_CONDITION_COLORS["treat"] if s in outlier_list else "#A8CADA"
              for s in features.index]
        sizes = [150 if s in outlier_list else 50 for s in features.index]

        for i, sample in enumerate(features.index):
            ax.scatter(embedding[i, 0], embedding[i, 1],
                      c=colors[i], s=sizes[i], alpha=PLOT_ALPHA,
                      edgecolors="black", linewidth=0.5)

            if colors[i] == "red":  # Annotate only outliers
                ax.annotate(sample, (embedding[i, 0], embedding[i, 1]),
                           fontsize=8, alpha=0.8, ha="center")

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_title(f"{outlier_type.capitalize()} Outliers (n={len(outlier_list)})")
        ax.legend(["Normal", "Outlier"], loc="best")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def plot_neighbors_network(embedding: np.ndarray, features: pd.DataFrame,
                          neighbors: Dict[str, List], outfile: Path,
                          axis_labels: Tuple[str, str] = ("UMAP 1", "UMAP 2")):
    """Plot UMAP with neighbor connections."""
    logger.info(f"Plotting neighbor networks -> {outfile}")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot edges for neighbors
    for sample, neighbor_list in neighbors.items():
        sample_idx = list(features.index).index(sample)
        for neighbor in neighbor_list[:3]:  # Show top 3 neighbors
            neighbor_idx = list(features.index).index(neighbor)
            ax.plot([embedding[sample_idx, 0], embedding[neighbor_idx, 0]],
                   [embedding[sample_idx, 1], embedding[neighbor_idx, 1]],
                   alpha=0.2, color="gray", linewidth=0.5)

    # Plot samples
    ax.scatter(embedding[:, 0], embedding[:, 1], s=100, alpha=PLOT_ALPHA,
              c="#5BAFD0", edgecolors="#1C4062", linewidth=0.5)

    # Annotate
    for i, sample in enumerate(features.index):
        ax.annotate(sample, (embedding[i, 0], embedding[i, 1]),
                   fontsize=8, alpha=0.8, ha="center")

    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.set_title("Sample Neighborhood Network (3-NN connections)")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def plot_global_metrics(embedding: np.ndarray, features: pd.DataFrame,
                       outfile: Path,
                       axis_labels: Tuple[str, str] = ("UMAP 1", "UMAP 2")):
    """Plot UMAP colored by various global metrics."""
    logger.info(f"Plotting global metrics -> {outfile}")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    aurora_map = LinearSegmentedColormap.from_list(
        "aurora", ["#A8CADA", "#5BAFD0", "#1C4062"]
    )

    # 1. Total modification rate (sum of all biotypes)
    if features.shape[1] > 0:
        total_mod = features.sum(axis=1)
        scatter = axes[0, 0].scatter(embedding[:, 0], embedding[:, 1],
                                     c=total_mod, cmap=aurora_map, s=100, alpha=PLOT_ALPHA)
        axes[0, 0].set_title("Total Modification Rate")
        plt.colorbar(scatter, ax=axes[0, 0])

    # 2. Feature diversity (number of non-zero features)
    feature_diversity = (features > 0).sum(axis=1)
    scatter = axes[0, 1].scatter(embedding[:, 0], embedding[:, 1],
                                c=feature_diversity, cmap=aurora_map, s=100, alpha=PLOT_ALPHA)
    axes[0, 1].set_title("Feature Diversity")
    plt.colorbar(scatter, ax=axes[0, 1])

    # 3. Local density (using 5-NN distance)
    neigh = NearestNeighbors(n_neighbors=5)
    neigh.fit(embedding)
    distances, _ = neigh.kneighbors(embedding)
    local_density = 1.0 / (distances[:, -1] + 1e-10)
    scatter = axes[1, 0].scatter(embedding[:, 0], embedding[:, 1],
                                c=local_density, cmap=aurora_map, s=100, alpha=PLOT_ALPHA)
    axes[1, 0].set_title("Local Density")
    plt.colorbar(scatter, ax=axes[1, 0])

    # 4. Maximum feature value
    max_feature = features.max(axis=1)
    scatter = axes[1, 1].scatter(embedding[:, 0], embedding[:, 1],
                                c=max_feature, cmap=aurora_map, s=100, alpha=PLOT_ALPHA)
    axes[1, 1].set_title("Max Feature Value")
    plt.colorbar(scatter, ax=axes[1, 1])

    # Common labels
    for ax in axes.flat:
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def save_results(embedding: np.ndarray, features: pd.DataFrame,
                cluster_labels: np.ndarray, outliers: Dict[str, List],
                neighbors: Dict[str, List], outdir: Path,
                metadata: Optional[pd.DataFrame] = None,
                axis_labels: Tuple[str, str] = ("UMAP1", "UMAP2")):
    """Save analysis results as TSV and JSON files."""
    logger.info(f"Saving results to {outdir}")

    outdir.mkdir(parents=True, exist_ok=True)

    # Save UMAP coordinates
    umap_df = pd.DataFrame(
        embedding,
        index=features.index,
        columns=list(axis_labels)
    )
    if metadata is not None:
        metadata = metadata.reindex(features.index)
        for column in ("condition", "cell_line"):
            if column in metadata:
                umap_df[column] = metadata[column].values
        if "condition" in metadata:
            umap_df["condition_label"] = metadata["condition"].map(CONDITION_LABELS).values
    umap_df.to_csv(outdir / "umap_coordinates.tsv", sep="\t")
    logger.info(f"Saved UMAP coordinates to umap_coordinates.tsv")

    membership_df = pd.DataFrame({
        "sample": features.index,
        "cluster_id": cluster_labels,
        "dbscan_membership": np.where(cluster_labels == -1, "noise", "cluster"),
    })
    membership_df.to_csv(outdir / "dbscan_membership.tsv", sep="\t", index=False)
    logger.info("Saved DBSCAN membership to dbscan_membership.tsv")

    # Save outlier classification
    outlier_classification = {}
    for outlier_type, outlier_list in outliers.items():
        for sample in features.index:
            if sample not in outlier_classification:
                outlier_classification[sample] = []
            if sample in outlier_list:
                outlier_classification[sample].append(outlier_type)

    outlier_df = pd.DataFrame([
        {"sample": s, "outlier_types": ",".join(outlier_classification.get(s, [""]))}
        for s in features.index
    ])
    outlier_df.to_csv(outdir / "outlier_classification.tsv", sep="\t", index=False)
    logger.info(f"Saved outlier classification to outlier_classification.tsv")

    # Save neighbors as JSON
    with open(outdir / "neighbors.json", "w") as f:
        json.dump(neighbors, f, indent=2)
    logger.info(f"Saved neighbor network to neighbors.json")

    # Save feature summary
    features.to_csv(outdir / "feature_matrix.tsv", sep="\t")
    logger.info(f"Saved feature matrix to feature_matrix.tsv")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="UMAP analysis of sRNA-TAPS modifications")
    parser.add_argument("--calls-dir", type=Path, default=None,
                       help="Directory with individual sample calls (e.g., 07.taps_calls)")
    parser.add_argument("--stringent-dir", type=Path, default=None,
                       help="Directory with stringent condition calls (e.g., 07e.stringent_calls)")
    parser.add_argument("--biotypes", default="rRNA,miRNA,tRNA,snoRNA,snRNA,piRNA,lncRNA,other",
                       help="Comma-separated list of biotypes")
    parser.add_argument("--samples-tsv", type=Path, default=None,
                       help="Optional sample metadata TSV")
    parser.add_argument("--outdir", required=True, type=Path,
                       help="Output directory for results and plots")
    parser.add_argument("--n-neighbors", type=int, default=15,
                       help="UMAP n_neighbors parameter")
    parser.add_argument("--eps", type=float, default=0.5,
                       help="DBSCAN eps parameter for clustering")
    parser.add_argument("--random-state", type=int, default=42,
                       help="Random seed")
    parser.add_argument("--method", choices=["umap", "pacmap", "pca"], default="umap",
                       help="Embedding method")
    parser.add_argument("--conditions", default=None,
                       help="Optional comma-separated internal conditions to include")
    parser.add_argument("--replicates-per-group", type=int, default=None,
                       help="Retain this many profiles closest to each condition-by-cell-line centroid")

    args = parser.parse_args()

    # Determine which calls directory to use
    calls_dir = None
    if args.calls_dir and args.calls_dir.exists():
        calls_dir = args.calls_dir
        logger.info(f"Using individual sample calls: {args.calls_dir}")
    elif args.stringent_dir and args.stringent_dir.exists():
        calls_dir = args.stringent_dir
        logger.info(f"Using stringent condition calls: {args.stringent_dir}")
    else:
        logger.error("Must provide either --calls-dir or --stringent-dir")
        return 1

    # Parse biotypes
    biotypes = [b.strip() for b in args.biotypes.split(",")]

    # Load and process data
    features = load_calls(calls_dir, biotypes, [])
    if features.empty:
        logger.error("Failed to load calls!")
        return 1

    features, metadata = engineer_features(features, calls_dir,
                                          biotypes, args.samples_tsv)
    if args.conditions and metadata is not None:
        selected_conditions = {condition.strip() for condition in args.conditions.split(",")}
        keep = metadata["condition"].isin(selected_conditions)
        features = features.loc[features.index.intersection(metadata.index[keep])]
        metadata = metadata.loc[features.index]
        logger.info("Restricted analysis to conditions: %s (%d samples)",
                    ", ".join(sorted(selected_conditions)), len(features))
        if features.empty:
            logger.error("No samples remain after applying --conditions")
            return 1
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.replicates_per_group is not None:
        features, metadata, selection_df = select_replicates_per_group(
            features, metadata, args.replicates_per_group
        )
        selection_df.to_csv(args.outdir / "replicate_selection.tsv", sep="\t", index=False)
        logger.info("Retained %d samples using %d replicates per condition/cell-line group",
                    len(features), args.replicates_per_group)
    run_statistics(features, metadata, args.outdir, random_state=args.random_state)

    # Compute UMAP
    embedding = compute_umap(features, n_neighbors=args.n_neighbors,
                            random_state=args.random_state, method=args.method)
    if args.method == "pca":
        scaled_features = StandardScaler().fit_transform(features)
        explained = PCA(n_components=2).fit(scaled_features).explained_variance_ratio_
        axis_labels = (f"PC1 ({explained[0] * 100:.1f}%)",
                       f"PC2 ({explained[1] * 100:.1f}%)")
        plot_title = "PCA of sample-level modification profiles"
    else:
        axis_labels = ("UMAP 1", "UMAP 2")
        plot_title = "Sample-level modification profiles"

    # Analyze structure
    clusters = detect_clusters(embedding, eps=args.eps)
    outliers = identify_outliers(embedding, features)
    neighbors = compute_local_neighbors(embedding, features)

    # Generate plots
    plot_umap_clusters(embedding, features, clusters,
                             args.outdir / "01_umap_clusters.png", metadata,
                             axis_labels, plot_title)
    plot_outliers(embedding, features, outliers,
                      args.outdir / "02_umap_outliers.png", axis_labels)
    plot_neighbors_network(embedding, features, neighbors,
                                  args.outdir / "03_umap_neighbors.png", axis_labels)
    plot_global_metrics(embedding, features,
                              args.outdir / "04_umap_metrics.png", axis_labels)

    # Save results
    save_results(embedding, features, clusters, outliers, neighbors, args.outdir,
                 metadata, tuple(label.replace(" ", "") for label in axis_labels))

    logger.info("UMAP analysis complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
