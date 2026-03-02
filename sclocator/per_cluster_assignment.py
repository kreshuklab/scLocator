import pandas as pd
from pathlib import Path
import argparse

def summarize_assignments_per_cluster(probability_file: str | Path):
    """
    Processes a cluster probability file to generate summaries for each cluster.

    This function creates two types of outputs:
    1. A folder ('per_cluster_assignments') containing a separate TSV file for each 
       cluster, listing all cells where that cluster was in the top 10.
    2. A single master TSV file ('master_top_10_ranks.tsv') showing the top-10
       ranks for all cells and all clusters in one table.
    """

    probability_file = Path(probability_file)
    
    try:
        df = pd.read_csv(probability_file, sep='\t', index_col='label_id')
        print(f"Successfully loaded '{probability_file}'.")
    except FileNotFoundError:
        print(f"Error: Could not find the file at '{probability_file}'.")
        print("Please make sure the path to the probability file is correct.")
        return

    most_probable_col = 'most probable cluster'
    if most_probable_col not in df.columns:
        print(f"Error: The file must contain a '{most_probable_col}' column.")
        return
        
    prob_matrix = df.drop(columns=[most_probable_col], errors='ignore')
    
    # Calculate the rank of each cluster for every cell
    ranks_df = prob_matrix.rank(axis=1, method='min', ascending=False)

    # Create the output directory
    output_dir = probability_file.parent / "per_cluster_assignments"
    output_dir.mkdir(exist_ok=True)
    print(f"Results will be saved in the '{output_dir}' directory.")

    # Generate Individual Per-Cluster Files
    processed_clusters = 0
    for cluster_name in prob_matrix.columns:
        cluster_ranks = ranks_df[cluster_name]
        top_10_mask = cluster_ranks <= 10
        
        if not top_10_mask.any():
            continue
            
        top_10_df = cluster_ranks[top_10_mask].reset_index()
        top_10_df.columns = ['label_id', 'rank']
        top_10_df['inverse_rank'] = 11 - top_10_df['rank']
        top_10_df['assignment_type'] = top_10_df['rank'].apply(lambda x: 'assigned' if x == 1 else 'top_10')
        top_10_df = top_10_df.sort_values(by='rank')
        
        safe_cluster_name = cluster_name.replace(' ', '_').replace('/', '_')
        output_filename = output_dir / f"{safe_cluster_name}_top_assignments.tsv"
        
        top_10_df.to_csv(output_filename, sep='\t', index=False)
        processed_clusters += 1

    print(f"Generated {processed_clusters} individual summary files.")
    
    # Create a copy of the ranks dataframe to modify
    master_ranks_df = ranks_df.copy()
    master_ranks_df[master_ranks_df > 10] = 0

    master_ranks_df = master_ranks_df.fillna(0)
    
    # Convert the entire dataframe to integers for a cleaner output
    master_ranks_df = master_ranks_df.astype(int)
    
    # Save the master table
    master_output_path = output_dir / "master_top_10_ranks.tsv"
    master_ranks_df.to_csv(master_output_path, sep='\t', index=True, index_label='label_id')
    
    print(f"Generated the master rank table at: '{master_output_path}'")
    print(f"\n All files have been created.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarizes cell assignments per cluster from a probability file."
    )
    
    parser.add_argument(
        "probability_file",
        nargs='?',
        default='cluster_probability.tsv',
        help="Path to the cluster probability TSV file. Defaults to 'cluster_probability.tsv' in the current directory."
    )
    
    args = parser.parse_args()
    
    summarize_assignments_per_cluster(args.probability_file)