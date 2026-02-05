import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np


def get_prob(expr: np.ndarray, profile: np.ndarray, gene_weights: np.ndarray) -> float:
    """
    Calculates the weighted log-likelihood of an expression vector given a profile.
    """
    positive_profile = profile
    negative_profile = 1 - profile
    positive_expr = expr
    negative_expr = 1 - expr

    # Calculate the individual likelihood term for each gene
    prob_vector = positive_expr * positive_profile + negative_expr * negative_profile

    # Switch to a weighted sum of log-probabilities
    # Add a small epsilon to prevent errors from log(0)
    log_likelihoods = np.log(prob_vector + 1e-10)

    # Apply the weights to each gene's log-likelihood
    weighted_log_likelihood = np.sum(log_likelihoods * gene_weights)

    # This function now returns a single log-likelihood score
    return weighted_log_likelihood


def expr_to_cluster_probability(
    expr: np.ndarray,
    cluster_prior_prob: float,
    cluster_profile: np.ndarray,
    gene_weights: np.ndarray,
) -> float:
    """
    Calculates the posterior probability of a cluster given an expression vector.
    """
    log_likelihood = get_prob(expr, cluster_profile, gene_weights)

    # Combine with the prior probability in log-space
    log_posterior = log_likelihood + np.log(cluster_prior_prob)

    # Return the final probability (not in log-space)
    return np.exp(log_posterior)


def get_individual_gene_probabilities(expr, profile):
    """
    Calculates the likelihood contribution (Term_g) for each gene
    without multiplying them together (i.e., P(E_g | C)).
    """
    positive_profile = profile
    negative_profile = 1 - profile

    positive_expr = expr
    negative_expr = 1 - expr

    # This is the vector of likelihoods (Term_g) for each gene
    prob_vector = positive_expr * positive_profile + negative_expr * negative_profile
    return prob_vector


def assign_clusters(
    spatial_expr_df,
    gene_list,
    profiles,
    output_dir,
    cluster_column,
    gene_weights: np.ndarray,
    run_diagnostics: bool = False,
    diagnostics_dir: str | Path = None,
):

    cluster_list = list(profiles.keys())
    result_df = pd.DataFrame(index=spatial_expr_df.index, columns=cluster_list)

    # Prior probability: each cluster has same prior probability
    cluster_prior_prob = 1 / len(cluster_list)

    for cluster_name in cluster_list:
        cluster_profile = profiles[cluster_name]

        bayesian_prob = [
            expr_to_cluster_probability(
                expr.to_numpy().flatten(),
                cluster_prior_prob,
                cluster_profile,
                gene_weights=gene_weights,
            )
            for _, expr in spatial_expr_df.iterrows()
        ]
        result_df[cluster_name] = bayesian_prob

    logging.info(f"Highest cluster probaility: {result_df.max()}")
    logging.info(f"Probability sum: {result_df.sum(axis=1)}")
    logging.info(f"Probability sum min: {result_df.sum(axis=1).min()}")

    # 1. IDENTIFY FAILURES AND GENERIC ASSIGNMENTS

    temp_assignment = result_df.idxmax(axis=1)
    sum_probs = result_df.sum(axis=1)

    nan_cells_mask = result_df.isna().any(axis=1)
    zero_sum_cells_mask = sum_probs <= 1e-300
    failed_cells_mask = nan_cells_mask | zero_sum_cells_mask
    failed_cell_ids = result_df[failed_cells_mask].index.tolist()

    if run_diagnostics:
        # --- Create the separate diagnostics directory ---
        diag_path = Path(diagnostics_dir)
        diag_path.mkdir(exist_ok=True)
        logging.info(
            f"Diagnostics enabled. Saving all diagnostic files to: '{diag_path}'"
        )

        if failed_cell_ids:
            logging.error(
                f"FATAL ASSIGNMENT FAILURE: Found {len(failed_cell_ids)} cells resulting in NaN or near-zero total probability."
            )

            failed_expressions = spatial_expr_df.loc[failed_cells_mask]
            failed_probabilities = result_df.loc[failed_cells_mask]
            failed_details_df = failed_expressions.join(
                failed_probabilities, rsuffix="_prob"
            )
            output_failed_path = diag_path / "failed_cell_details.tsv"
            failed_details_df.to_csv(
                output_failed_path, sep="\t", index=True, index_label="label_id"
            )
            logging.info(
                f"Saved detailed diagnostics for {len(failed_details_df)} failed cells to '{output_failed_path}'."
            )

        # Save diagnostics for generic clusters
        generic_cells_mask = temp_assignment.isin(["zero", "autofluorescence"])
        df_generic_unnormalized = result_df[
            generic_cells_mask & ~failed_cells_mask
        ].copy()

        if not df_generic_unnormalized.empty:
            df_generic_unnormalized["assigned_cluster"] = temp_assignment[
                df_generic_unnormalized.index
            ]
            output_diagnostic_path = diag_path / "generic_cluster_diagnostics.tsv"
            df_generic_unnormalized.to_csv(
                output_diagnostic_path, sep="\t", index=True, index_label="label_id"
            )
            logging.info(
                f"Saved diagnostics for {len(df_generic_unnormalized)} generic cells."
            )

        # 2. CREATE MASTER DIAGNOSTIC SAMPLE LIST
        MAX_SAMPLES_PER_CLUSTER = 5
        diagnostic_cell_ids = list(set(failed_cell_ids))
        successful_cells_assignment = temp_assignment[~failed_cells_mask]

        if not successful_cells_assignment.empty:
            ASSIGNMENT_COL_TEMP = "assigned_cluster_id"
            temp_df = successful_cells_assignment.to_frame(
                name=ASSIGNMENT_COL_TEMP
            ).reset_index()
            sampled_successful_ids = (
                temp_df.groupby(ASSIGNMENT_COL_TEMP)["label_id"]
                .apply(lambda x: x.sample(min(len(x), MAX_SAMPLES_PER_CLUSTER)))
                .tolist()
            )
            diagnostic_cell_ids = list(
                set(diagnostic_cell_ids + sampled_successful_ids)
            )

        diagnostic_cell_ids = diagnostic_cell_ids[
            : MAX_SAMPLES_PER_CLUSTER * len(cluster_list)
        ]

        # 3. SAVE COMPREHENSIVE GENE LIKELIHOOD DIAGNOSTICS

        if diagnostic_cell_ids:
            df_diagnostics_gene_likelihood = []

            for cell_id in diagnostic_cell_ids:
                cell_assignment = temp_assignment.loc[cell_id]
                cell_expr = spatial_expr_df.loc[cell_id].to_numpy()

                # Calculate the likelihood vector for this cell against *ALL* clusters
                for cluster_name in cluster_list:
                    profile = profiles[cluster_name]
                    gene_likelihoods = get_individual_gene_probabilities(
                        cell_expr, profile
                    )

                    df_row = pd.DataFrame(
                        {
                            "label_id": cell_id,
                            "assigned_cluster": cell_assignment,
                            "tested_cluster": cluster_name,
                            "gene_name": gene_list,
                            "likelihood_term": gene_likelihoods,
                        }
                    )
                    df_diagnostics_gene_likelihood.append(df_row)

            if df_diagnostics_gene_likelihood:
                df_diagnostics_final_gene = pd.concat(
                    df_diagnostics_gene_likelihood, ignore_index=True
                )

                if cluster_column == "curated_high_res_annotations":
                    logging.info("Splitting likelihood diagnostics by tested cluster.")
                    output_dir_likelihood = diag_path / "gene_likelihood_vectors"
                    output_dir_likelihood.mkdir(parents=True, exist_ok=True)

                    for cluster_name, df_group in df_diagnostics_final_gene.groupby(
                        "tested_cluster"
                    ):
                        safe_cluster_name = cluster_name.replace(" ", "_").replace(
                            "/", "_"
                        )
                        output_path = (
                            output_dir_likelihood
                            / f"likelihood_{safe_cluster_name}.tsv"
                        )
                        df_group.to_csv(output_path, sep="\t", index=False)

                else:
                    output_likelihood_path = (
                        diag_path / "all_cluster_gene_likelihood_diagnostic.tsv"
                    )
                    df_diagnostics_final_gene.to_csv(
                        output_likelihood_path, sep="\t", index=False
                    )
                    logging.info(
                        f"Saved comprehensive gene likelihood diagnostics to: {output_likelihood_path}"
                    )

    # Normalize the probabilities to sum to 1
    result_df = result_df.div(result_df.sum(axis=1), axis=0)

    logging.info(f"Highest cluster probaility: {result_df.max()}")
    logging.info(f"Probability sum min: {result_df.sum(axis=1).min()}")

    result_df["most probable cluster"] = result_df.idxmax(axis=1)

    result_df.to_csv(output_dir / "cluster_probability.tsv", sep="\t")
