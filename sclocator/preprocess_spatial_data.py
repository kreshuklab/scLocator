import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import h5py
import logging

from pathlib import Path
from typing import List


def percentile_norm(expr, pmin, pmax, eps=1e-10):
    """Percentile normalization of an expression vector. Values are clipped to [0, 1] range
    Args:
        expr: expression vector
        pmin: in percents
        pmax: in percents
        eps: _description_. Defaults to 1e-10.

    Returns:
        _description_
    """

    pmin = np.percentile(expr, pmin)
    pmax = np.percentile(expr, pmax)
    return np.clip(
        ((expr.astype(np.float32) - pmin) / (pmax - pmin + eps)), a_min=0, a_max=1
    )


def plot_cumulative_expression(expr, gene_name):
    expr_sorted = np.sort(expr)
    # plt.figure()
    plt.plot(expr_sorted, label=gene_name)
    plt.legend(bbox_to_anchor=(1.1, 1.05), ncol=3)
    plt.title("Cumulative normalized distribution")


def preprocess_spatial_expr_df(
    output_dir: Path,
    name: str,
    df_path: Path,
    bin_thresh: float,
    visualize: bool,
    apply_norm: bool,
    pmin: float,
    pmax: float,
    metadata_column: str,
):
    df_expr = pd.read_csv(df_path, sep="\t")

    assert (
        "label_id" in df_expr.columns
    ), "Spatial expression table should contain 'label_id' column"

    df_expr = df_expr.set_index("label_id")
    # Group the duplicate columns (with OR)
    if df_expr.columns.duplicated().any():
        df_expr = df_expr.groupby(df_expr.columns, axis=1).max()

    logging.info(f"Quality control and preprocessing of the spatial data {name}")

    logging.info("\n" + repr(df_expr.describe().T))

    # Percentile norm
    df_norm_expr = df_expr.copy(deep=True)
    df_norm_expr = df_norm_expr.apply(lambda x: percentile_norm(x, pmin, pmax), axis=0)

    # Binarize
    if apply_norm:
        df_bin_expr = df_norm_expr.copy(deep=True)
        df_bin_expr[df_norm_expr > bin_thresh] = 1
        df_bin_expr[df_norm_expr <= bin_thresh] = 0
    else:
        df_bin_expr = df_expr.copy(deep=True)
        df_bin_expr[df_expr >= bin_thresh] = 1
        df_bin_expr[df_expr < bin_thresh] = 0

    logging.info("Normalized expression df stats\n" + repr(df_norm_expr.describe().T))
    logging.info("Binarized expression df stats\n" + repr(df_bin_expr.describe().T))

    if False:
        plt.figure(figsize=(10, 1 * len(df_expr.columns)))
        sns.stripplot(df_expr, jitter=0.4, alpha=0.6, orient="h")
        plt.title("Distribution of raw counts")
        plt.savefig(output_dir / "raw_count_stripplot.png", dpi=300)
        plt.close("all")

        plt.figure(figsize=(10, 1 * len(df_expr.columns)))
        sns.violinplot(df_expr, alpha=0.6, orient="h", cut=0)
        plt.title("Distribution of raw counts")
        plt.savefig(output_dir / "raw_count_violinplot.png", dpi=300)
        plt.close("all")

        # Normalized expression per gene
        cumulative_distr_dir = output_dir / "cumulative_expr"
        cumulative_distr_dir.mkdir(exist_ok=True)

        for col in df_norm_expr.columns:
            plt.figure(figsize=(10, 10))
            plot_cumulative_expression(df_norm_expr[col], col)
            plt.tight_layout()
            plt.savefig(
                cumulative_distr_dir / f"norm_count_cumulative_{col}.png", dpi=300
            )
            plt.close("all")

        # Nonzero cells per gene
        plt.figure(figsize=(1 * len(df_expr.columns), 5))
        plt.bar(df_bin_expr.columns, df_bin_expr.sum(axis=0))
        plt.title("Nonzero cells per gene in binarized data")
        plt.xlabel("genes")
        plt.ylabel("Number of cells")
        plt.savefig(output_dir / "bin_nonzero_per_gene.png", dpi=300)
        plt.close("all")

        # Nonzero genes per cell
        plt.figure(figsize=(0.5 * len(df_expr.columns), 5))
        plt.bar(
            df_bin_expr.sum(axis=1).value_counts().index,
            df_bin_expr.sum(axis=1).value_counts(),
        )
        plt.xlabel("Number of nonzero genes")
        plt.ylabel("Number of cells")
        plt.savefig(output_dir / "bin_nonzero_per_cell.png", dpi=300)
        plt.close("all")

    return df_bin_expr


def harmonize_gene_names(metadata, expr_df, atlas_col, spatial_col):
    rename_dict = {
        spatial: atlas
        for spatial, atlas in zip(metadata[spatial_col], metadata[atlas_col])
    }

    missing_genes = [
        col
        for col in expr_df.columns
        if not (col in list(metadata[spatial_col]))
    ]
    print(list(metadata[spatial_col]))
    for col in expr_df.columns:
        logging.info(col in list(metadata[spatial_col]))
    logging.warning(
        f"Genes in spatial source without mapping to scRNAseq atlas: {missing_genes}"
    )
    columns = [
        col for col in expr_df.columns if col in list(metadata[spatial_col])
    ]

    renamed_expr_df = expr_df[columns].copy()
    renamed_expr_df = renamed_expr_df.rename(columns=rename_dict)

    return renamed_expr_df


def join_spatial_expr_dfs(
    dfs: List, df_names: List, metadata: pd.DataFrame, atlas_col, gene_name_col: str
):
    # Check that index is the same
    logging.info("\n\n\nJoin spatial dataframes")
    logging.info(f"name\t\tN genes\t\tN cells")
    for df, name in zip(dfs, df_names):
        logging.info(f"{name}\t{len(df.columns)}\t{len(df)}")
    if len(dfs) == 1:
        return dfs[0]
    else:
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = pd.concat([merged_df, df], axis=1, join="inner")

    # Group the duplicate columns (with OR)
    gene_name_dict = {
        spatial: atlas
        for spatial, atlas in zip(metadata[atlas_col], metadata[gene_name_col])
    }
    if merged_df.columns.duplicated().any():
        duplicated_columns = merged_df.columns[merged_df.columns.duplicated()]
        logging.info(f"Duplicated columns")
        for col in duplicated_columns:
            logging.info(f"{col} \t {gene_name_dict[col]}")
        merged_df = merged_df.T.groupby(merged_df.columns).max().T

    logging.info(f"Merged spatial dataframe")
    logging.info(f"name\t\tN genes\t\tN cells")
    logging.info(f"merged\t{len(merged_df.columns)}\t{len(merged_df)}")

    if np.any(merged_df.isna()):
        logging.warning(f"Merged spatial dataframe has NaNs")
        merged_df.fillna(0, inplace=True)

    return merged_df


def calculate_spatial_specificity_weights(
    spatial_expr_df: pd.DataFrame, min_weight: float = 0.25, max_weight: float = 1.75
) -> np.ndarray:
    """
    Calculates a weight for each gene based on how many nuclei express it.
    """
    logging.info("Calculating spatial gene specificity weights...")

    num_genes = spatial_expr_df.shape[1]
    num_nuclei = spatial_expr_df.shape[0]
    nuclei_counts = spatial_expr_df.sum(axis=0).to_numpy()
    ubiquity_score = nuclei_counts / num_nuclei
    specificity_score = 1 - ubiquity_score
    spatial_weights = min_weight + specificity_score * (max_weight - min_weight)

    logging.info(
        f"Spatial weights calculated, ranging from {np.min(spatial_weights):.2f} to {np.max(spatial_weights):.2f}"
    )

    return spatial_weights


def calculate_gene_specificity_weights(
    profiles: dict,
    gene_list: list,
    expression_threshold: float = 0.1,
    min_weight: float = 0.25,
    max_weight: float = 1.75,
) -> np.ndarray:
    """
    Calculates a weight for each gene based on its expression specificity across clusters.
    ... (rest of the docstring) ...
    """
    logging.info("Calculating gene specificity weights...")

    num_genes = len(gene_list)
    num_clusters = len(profiles)
    profile_matrix = np.array(
        list(profiles.values())
    )  # Shape: (num_clusters, num_genes)
    is_expressed_matrix = profile_matrix > expression_threshold
    gene_expression_counts = is_expressed_matrix.sum(axis=0)

    # Calculate "ubiquity score" (0 = perfectly specific, 1 = in all clusters)
    # We use max(1, ...) to avoid division by zero if there's only one cluster
    ubiquity_score = (gene_expression_counts - 1) / max(1, num_clusters - 1)
    ubiquity_score = np.clip(ubiquity_score, 0, 1)  # Ensure score is between 0 and 1

    # Invert to get "specificity score" (1 = perfectly specific, 0 = ubiquitous)
    specificity_score = 1 - ubiquity_score

    # Map specificity score to the desired weight range
    gene_weights = min_weight + specificity_score * (max_weight - min_weight)

    logging.info(
        f"Gene weights calculated, ranging from {np.min(gene_weights):.2f} to {np.max(gene_weights):.2f}"
    )

    return gene_weights, gene_expression_counts