import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np
import pickle


def get_atlas(
    adata_path,
    cluster_column,
    normalize_expression,
    layer,
    metadata_column,
    bin_thresh,
    output_dir,
    visualize,
    gene_list,
    metadata,
    vis_name_metadata_column,
):
    adata = sc.read_h5ad(adata_path)
    logging.info(adata)


    if layer != "X":
        adata.layers["counts"] = adata.layers.pop(layer)
        adata.X = adata.layers["counts"].copy()

    # Normalize
    if normalize_expression:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if visualize:
        sc.pl.umap(adata, color=cluster_column, size=2, show=False)
        plt.savefig(output_dir / "umap.png", dpi=300, bbox_inches="tight")
        plt.close("all")

    logging.info(
        "Number of cells per cluster/n" + str(adata.obs[cluster_column].value_counts())
    )
    # Assign metadata
    metadata = metadata[metadata[metadata_column].isin(gene_list)]
    genes_lookup = dict(
        zip(metadata[metadata_column], metadata[vis_name_metadata_column])
    )


    missing_genes = set(gene_list).difference(set(adata.var_names.to_list()))
    logging.info("Genes present in metadata but missing from the atlas")
    for gene in missing_genes:
        logging.info(f"Missing gene: {gene}")
    assert not set(gene_list).difference(
        set(adata.var_names.to_list())
    ), "Some genes are not present in the metadata table"

    
    subset_adata = adata[:, adata.var.index.isin(gene_list)].copy()
    subset_adata.var[metadata_column] = subset_adata.var_names
    subset_adata.var["vis_name"] = [
        genes_lookup[name] for name in subset_adata.var_names
    ]

    if bin_thresh == "median":
        thresholds = []
        for name in subset_adata.var.index:
            expr = subset_adata[:, name].X.toarray().T[0]
            expr_sorted = np.sort(expr)
            # Calculate the median of non-zero expression values
            thr = np.median(expr_sorted[expr_sorted > 0])
            thresholds.append(thr)
    else:
        thresholds = [bin_thresh] * len(subset_adata.var_names)

    if visualize:
        logging.info("Generating cumulative expression plots...")
        for i, name in enumerate(subset_adata.var.index):
            print(f"Plot cumulative expression {name} | {genes_lookup[name]}")
            plt.figure(figsize=(6, 5))
            expr = subset_adata[:, name].X.toarray().T[0]
            expr_sorted = np.sort(expr)

            # Use the pre-calculated threshold for the plot
            thr = thresholds[i]

            plt.plot(expr_sorted, label=name)
            plt.hlines(y=thr, xmin=0, xmax=len(expr_sorted), color="red")
            plt.legend()
            plt.title(f"Cumulative distribution of expression levels {name} | {genes_lookup[name]}")
            plt.ylabel("Normalized expression level")
            plt.xlabel("Cell number")
            plt.savefig(
                output_dir
                / f"{subset_adata.var.loc[name, 'vis_name']}_norm_count_cumulative.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close("all")

    logging.info(subset_adata)

    subset_adata.write(output_dir / "scrnaseq_norm_expr.h5ad")

    subset_adata_bin = subset_adata.copy()
    subset_adata_bin.X = (subset_adata_bin.X > np.array(thresholds)).astype(np.int16)
    subset_adata_bin.write(output_dir / "scrnaseq_bin_expr.h5ad")

    return subset_adata_bin



def get_cluster_profiles(adata, cluster_column, output_dir, norm_per_gene, base_prob: float = 0.001, delta = 0.0000001):

    adata_aggregated = sc.get.aggregate(adata, by=cluster_column, func="count_nonzero")
    logging.info(adata_aggregated)

    cluster_sizes = adata.obs[cluster_column].value_counts()
    cluster_names = adata.obs[cluster_column].unique()
    logging.info(cluster_sizes)
    profiles = {}
    max_expr_prob = np.zeros(adata.shape[1])
    min_expr_prob = np.ones(adata.shape[1])
    for cluster_name in cluster_names:
        nonzero_counts = adata_aggregated[cluster_name, :].layers["count_nonzero"]

        # This ensures P is never exactly 0 or 1.
        profiles[cluster_name] = (nonzero_counts + base_prob) / (
            cluster_sizes[cluster_name] + 2 * base_prob
        )
        profiles[cluster_name] = np.array(profiles[cluster_name]).flatten()

        max_expr_prob = np.maximum(max_expr_prob, profiles[cluster_name])
        min_expr_prob = np.minimum(min_expr_prob, profiles[cluster_name])



    if norm_per_gene:
        norm_profiles = {"min": min_expr_prob, "max": max_expr_prob}
        output_dir_minmax = output_dir / "profile_normalization"
        output_dir_minmax.mkdir(parents=True, exist_ok=True)
        visualize_profiles(list(adata.var.vis_name), norm_profiles, output_dir_minmax)
        for key in profiles:
            # profiles[key] = (profiles[key] - min_expr_prob) * (1 - 2 * base_prob) / (max_expr_prob - min_expr_prob) + base_prob
            profiles[key] = (profiles[key] - min_expr_prob) * (1 - 0.1 - base_prob) / (max_expr_prob - min_expr_prob) + base_prob
    # Create empty profile for the nuclei where nothing was measured
    profiles["zero"] = np.zeros(adata.shape[1]) + delta

    # Create profile that expresses everything for the autofluorescent areas
    profiles["autofluorescence"] = np.ones(adata.shape[1]) - 0.01

    output_file = output_dir / "profiles_dictionary.pkl"
    try:
        with open(output_file, "wb") as f:
            pickle.dump(profiles, f)
        logging.info(
            f"Successfully saved cluster profiles dictionary to: {output_file}"
        )
    except Exception as e:
        logging.error(f"Failed to save cluster profiles via pickle: {e}")

    return profiles


def visualize_profiles(gene_list, profiles, output_dir):
    """
    Generates and saves the gene expression profile for each cluster as a separate PNG file.

    Args:
        gene_list (list): List of gene names (for the x-axis).
        profiles (dict): Dictionary where keys are cluster names and values are
                         the expression profile arrays (0-1).
        output_dir (Path): The directory path (e.g., log_dir / "profiles") to save the plots.
    """
    logging.info(
        f"Generating and saving individual expression profile plots for {len(profiles)} clusters..."
    )

    # Ensure the output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Iterate through each cluster profile
    for cluster_name, profile_data in profiles.items():

        # 1. Create a new figure for each cluster
        # Use a reasonable size that accommodates all genes on the x-axis
        # Size suggestion: 0.15 inches per gene + margin, but we'll use a fixed wide size.
        fig, ax = plt.subplots(figsize=(15, 6))

        # 2. Plot the profile data
        ax.bar(gene_list, profile_data, color="darkblue", alpha=0.9)

        # 3. Set limits and labels
        ax.set_ylim(0, 1)
        ax.set_title(f"Cluster Profile: {cluster_name}", fontsize=14)
        ax.set_ylabel("Proportion of Cells Expressing Gene")
        ax.set_xlabel("Gene")

        # Rotate x-axis ticks for readability (important for many genes)
        ax.tick_params("x", rotation=90, labelsize=6)

        # 4. Save the figure
        # Sanitize the filename to handle spaces/slashes in cluster names
        safe_name = cluster_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
        plot_path = output_dir / f"profile_{safe_name}.png"

        # Use tight_layout only for this specific figure
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close(fig)  # Close the figure to free memory

    logging.info(
        f"Successfully saved {len(profiles)} individual cluster profile plots to {output_dir}"
    )

    output_dir = output_dir / "prob_per_gene"
    output_dir.mkdir(parents=True, exist_ok=True)

    profile_names = profiles.keys()
    for idx, gene in enumerate(gene_list):

        gene_prob_profile = np.zeros(len(profile_names))
        for j, name in enumerate(profile_names):
            gene_prob_profile[j] = profiles[name][idx]
        # 1. Create a new figure for each cluster
        # Use a reasonable size that accommodates all genes on the x-axis
        # Size suggestion: 0.15 inches per gene + margin, but we'll use a fixed wide size.
        fig, ax = plt.subplots(figsize=(15, 6))

        # 2. Plot the profile data
        ax.bar(profile_names, gene_prob_profile, color="darkblue", alpha=0.9)

        # 3. Set limits and labels
        ax.set_ylim(0, 1)
        ax.set_title(f"Gene name: {gene}", fontsize=14)
        ax.set_ylabel("Proportion of Cells Expressing Gene")
        ax.set_xlabel("Gene")

        # Rotate x-axis ticks for readability (important for many genes)
        ax.tick_params("x", rotation=90, labelsize=6)

        # 4. Save the figure
        # Sanitize the filename to handle spaces/slashes in cluster names
        plot_path = output_dir / f"profile_{gene}.png"

        # Use tight_layout only for this specific figure
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close(fig)  # Close the figure to free memory