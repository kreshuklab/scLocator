import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import sys
import pickle
import click
from pathlib import Path

from sclocator.bayesian_prob import get_individual_gene_probabilities

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")



# --- HELPER FUNCTIONS ---

def load_data(base_dir: Path, target_id: float):
    """Loads necessary TSVs and Profiles from the data repository."""
    expr_path = base_dir / "spatial_expr_df_binarized.tsv"
    assign_path = base_dir / "cluster_prob" / "cluster_probability.tsv"
    profile_path = base_dir / "expression_profiles" / "profiles_dictionary.pkl"
    
    # Load expression
    df_expr = pd.read_csv(expr_path, sep='\t', index_col='label_id')
    df_expr.index = df_expr.index.astype(float)
    
    # Load assignments
    df_assign = pd.read_csv(assign_path, sep=r'\s+', index_col='label_id')
    df_assign.index = df_assign.index.astype(float)
    
    # Load profiles
    with open(profile_path, 'rb') as f:
        profiles = pickle.load(f)
        
    return df_expr, df_assign, profiles

def plot_diagnostic_heatmap_OLD(df, cell_id, assigned, out_path, top_n_genes=40):
    """Generates the heatmap from calculated likelihoods."""
    # Pivot logic
    df['log_likelihood'] = np.log10(df['likelihood_term'] + 1e-10) 
    pivot_df = df.pivot(index='gene_name', columns='tested_cluster', values='log_likelihood')
    
    # Filter for discriminatory genes
    gene_variance = pivot_df.std(axis=1).sort_values(ascending=False)
    top_genes = gene_variance.head(top_n_genes).index
    plot_data = pivot_df.loc[top_genes]
    
    plt.figure(figsize=(14, 12))
    sns.heatmap(plot_data, cmap='RdYlGn', center=-2, 
                cbar_kws={'label': 'Log10 Likelihood (Green=Match, Red=Mismatch)'})
    
    plt.title(f"Diagnostic for Cell {cell_id}\nTop {top_n_genes} Discriminatory Genes", fontsize=16)
    plt.xticks(rotation=45, ha='right')
    
    # Highlight assigned cluster in bold
    for label in plt.gca().get_xticklabels():
        if label.get_text() == str(assigned):
            label.set_weight('bold')
            label.set_color('darkgreen')
            
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    logging.info(f"Heatmap saved to: {out_path}")


def plot_diagnostic_heatmap(df, cell_id, assigned, out_path, compare_targets=None, top_n_genes=25):
    """
    Generates a three-panel diagnostic heatmap:
    1. SUPPORTING: High likelihood in the winner.
    2. DISCRIMINATORY: High variance across all tested clusters.
    3. MISMATCH: High Likelihood Gap (Winner Match vs Competitor Mismatch).
    """
    # 1. Prepare Data

    raw_pivot_df = df.pivot(index='gene_name', columns='tested_cluster', values='likelihood_term')

    df['log_likelihood'] = np.log10(df['likelihood_term'] + 1e-10) 
    pivot_df = df.pivot(index='gene_name', columns='tested_cluster', values='log_likelihood')
    
    # Ensure assigned cluster is the baseline (first column)
    competitors = [c for c in pivot_df.columns if c != assigned]
    pivot_df = pivot_df[[assigned] + competitors]
    raw_pivot_df = raw_pivot_df[[assigned] + competitors] # <-- Keep raw aligned with the colors!

    # 2. Select Genes for each panel
    # Panel 1: Top likelihood in winner
    supporting_genes = pivot_df[assigned].sort_values(ascending=False).head(top_n_genes).index
    
    # Panel 2: Highest variance (Discriminatory)
    gene_variance = pivot_df.std(axis=1).sort_values(ascending=False)
    discriminatory_genes = gene_variance.head(top_n_genes).index
    
    # Panel 3: highest mismatch with comparison cluster

    # Ensure compare_targets is a list
    if compare_targets is None:
        compare_targets = []

    # Decide who Panel 3 should focus on
    if len(compare_targets) > 0:
        # Focus ONLY on the clusters you explicitly typed in the terminal
        panel3_targets = [c for c in compare_targets if c in pivot_df.columns and c != assigned]
    else:
        # Fallback to the Top 10 if you didn't use --compare_to
        panel3_targets = competitors

    # Panel 3: Absolute Mismatch (The Surgical Post-Mortem)
    if len(panel3_targets) > 0:
        mismatch_genes = pivot_df[panel3_targets].min(axis=1).sort_values().head(top_n_genes).index
    else:
        mismatch_genes = discriminatory_genes

    # 3. Create Triple-Panel Plot
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(32, 12))

    # Formatting helper for heatmaps
    heatmap_kwargs = {
        'cmap': 'RdYlGn', 
        'vmin': -2.0,  
        'vmax': 0.0,   
        'cbar_kws': {'label': 'Log10 Likelihood (Mismatch strength)', 'orientation': 'horizontal'}
    }

    # Plot Panels
    # sns.heatmap(pivot_df.loc[supporting_genes], ax=ax1, **heatmap_kwargs)
    sns.heatmap(pivot_df.loc[supporting_genes], ax=ax1, 
                annot=raw_pivot_df.loc[supporting_genes], fmt=".2f", annot_kws={"size": 8}, 
                **heatmap_kwargs)
    ax1.set_title(f"1. Genes Supporting Matches for {assigned}", fontsize=14, fontweight='bold')

    # sns.heatmap(pivot_df.loc[discriminatory_genes], ax=ax2, **heatmap_kwargs)
    sns.heatmap(pivot_df.loc[discriminatory_genes], ax=ax2, 
                annot=raw_pivot_df.loc[discriminatory_genes], fmt=".2f", annot_kws={"size": 8}, 
                **heatmap_kwargs)
    ax2.set_title("2. Most Variable Genes", fontsize=14, fontweight='bold')

    # sns.heatmap(pivot_df.loc[mismatch_genes], ax=ax3, **heatmap_kwargs)
    sns.heatmap(pivot_df.loc[mismatch_genes], ax=ax3, 
                annot=raw_pivot_df.loc[mismatch_genes], fmt=".2f", annot_kws={"size": 8}, 
                **heatmap_kwargs)
    ax3.set_title("3. Biggest Mismatches for Competitor Clusters", fontsize=14, fontweight='bold')

    # sns.heatmap(..., annot=raw_pivot_df.loc[supporting_genes], fmt=".2f")

    # X-Axis and Highlighting logic
    for ax in [ax1, ax2, ax3]:
        plt.sca(ax)
        plt.xticks(rotation=45, ha='right')
        
        # Highlight assigned cluster in bold
        for label in ax.get_xticklabels():
            if label.get_text() == str(assigned):
                label.set_weight('bold')
                label.set_color('darkgreen')
        
        ax.set_ylabel("")

    fig.suptitle(f"Triple-Panel Inspection: Cell/Nucleus {cell_id}", fontsize=24, y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plt.savefig(out_path, dpi=300)
    logging.info(f"Triple-panel plot saved to: {out_path}")


# --- MAIN CLICK COMMAND ---

@click.command()
@click.argument("nucleus_id", type=float)
@click.argument("base_dir", type=click.Path(exists=True))
@click.option("--top_n", default=10, help="Include top N probable clusters")
@click.option("--compare_to", multiple=True, help="Specific clusters to compare")
@click.option("--gene_id_col", default="XLOC", help="Column name for raw gene IDs")
@click.option("--gene_name_col", default="Platybrowser name", help="Column name for readable gene names")
def main(nucleus_id, base_dir, top_n, compare_to, gene_id_col, gene_name_col):
    """Unified script to calculate likelihoods and plot diagnostic heatmaps."""
    BASE_DIR = Path(base_dir)
    
    # 1. Load Data
    df_smfish, df_assignments, profiles = load_data(BASE_DIR, nucleus_id)
    
    # 2. Determine clusters to test
    probs_row = df_assignments.loc[[nucleus_id]].drop(columns=['most probable cluster'], errors='ignore').transpose()
    probs_row.columns = ['prob']
    probs_row['prob'] = pd.to_numeric(probs_row['prob'], errors='coerce')
    probs_row = probs_row.dropna(subset=['prob'])
    
    # 2. Derive the top N clusters (highest probability first)
    sorted_probs = probs_row.sort_values(by='prob', ascending=False)
    top_clusters = sorted_probs.head(top_n).index.tolist()

    assigned_cluster = sorted_probs.index[0]

    final_target_list = list(set(top_clusters + list(compare_to) + [assigned_cluster]))

    # 2.5 Gene name mapping
    raw_gene_list = df_smfish.columns.tolist()
    METADATA_PATH = Path("data/scrnaseq_atlas/genes_metadata.tsv")
    
    if METADATA_PATH.exists():
        logging.info(f"Loading gene metadata from {METADATA_PATH}...")
        meta_df = pd.read_csv(METADATA_PATH, sep='\t')
        
        # Create dictionary: { 'XLOC_ID': 'Gene_Symbol' }
        # Adjust 'xloc' and 'gene_symbol' below if your column headers are different
        # name_map = dict(zip(meta_df['XLOC'], meta_df['Platybrowser name']))
        name_map = dict(zip(meta_df[gene_id_col], meta_df[gene_name_col]))
        
        # Translate the list, keeping original if no match found
        gene_list = [name_map.get(g, g) for g in raw_gene_list]
        logging.info("Gene names successfully mapped.")
    else:
        logging.warning(f"Metadata not found at {METADATA_PATH}. Using raw column names.")
        gene_list = raw_gene_list
    
    # 3. Calculate Likelihoods
    cell_expr = df_smfish.loc[nucleus_id].to_numpy()
    # gene_list = df_smfish.columns.tolist()                -----!!!!!! MAYBE CHANGE BACK LATER???
    results = []

    for cluster_name in final_target_list:
        if cluster_name in profiles:
            likes = get_individual_gene_probabilities(cell_expr, profiles[cluster_name])
            results.append(pd.DataFrame({
                'label_id': nucleus_id,
                'assigned_cluster': assigned_cluster,
                'tested_cluster': cluster_name,
                'gene_name': gene_list,
                'likelihood_term': likes
            }))

    df_final = pd.concat(results, ignore_index=True)
    
    # 4. Save TSV and Plot
    safe_id = str(nucleus_id).replace('.', '_')
    tsv_out = BASE_DIR / "targeted_inspections" / f"cell_{safe_id}_vs_targets.tsv"
    plot_out = BASE_DIR / "targeted_inspections" / f"cell_{safe_id}_heatmap.png"
    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    
    df_final.to_csv(tsv_out, sep='\t', index=False)
    plot_diagnostic_heatmap(df_final, nucleus_id, assigned_cluster, plot_out, compare_targets=list(compare_to))

if __name__ == "__main__":
    main()