import pandas as pd
import numpy as np
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import click
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# --- HELPER FUNCTIONS ---
def generate_confidence_summary(base_dir: Path):
    prob_path = base_dir / "cluster_prob" / "cluster_probability.tsv"
    if not prob_path.exists(): return None
    
    df = pd.read_csv(prob_path, sep='\t', index_col='label_id')
    probs_only = df.drop(columns=['most probable cluster', 'cluster'], errors='ignore')
    max_probs = probs_only.max(axis=1)

    margins = probs_only.apply(lambda row: row.nlargest(2).values[0] - row.nlargest(2).values[1] if len(row.nlargest(2)) >= 2 else 0.0, axis=1)

    summary_df = pd.DataFrame({
        'label_id': df.index,
        'assigned_cluster': df['most probable cluster'] if 'most probable cluster' in df.columns else probs_only.idxmax(axis=1),
        'max_probability': max_probs,
        'confidence_margin': margins
    })

    summary_df.to_csv(base_dir / "cell_confidence_summary.tsv", sep='\t', index=False)
    return summary_df

# def plot_cluster_confidence(summary_df, base_dir: Path):
    stats = summary_df.groupby('assigned_cluster').agg({'max_probability': ['mean', 'min', 'max'], 'confidence_margin': ['mean', 'min', 'max']}).sort_index() 
    clusters = stats.index
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(20, len(clusters) * 0.25), 16), sharex=True)

    def plot_metric(ax, metric, color, ecolor, title, ylabel):
        mean_val = stats[(metric, 'mean')]
        ax.errorbar(clusters, mean_val, yerr=[mean_val - stats[(metric, 'min')], stats[(metric, 'max')] - mean_val], fmt='o', color=color, ecolor=ecolor, elinewidth=1, capsize=2, markersize=4)
        ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
        ax.set_title(title, fontsize=18)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(axis='both', linestyle=':', alpha=0.5)

    plot_metric(ax1, 'max_probability', 'navy', 'skyblue', 
                f'Assignment Confidence Audit ({len(clusters)} Cell Types) - Mean ± Std Dev', 'Mean Best-Match Score')
    plot_metric(ax2, 'confidence_margin', 'darkred', 'salmon', '', 'Confidence Margin (1st vs 2nd)')

    plt.xticks(rotation=90, fontsize=9)
    plt.xlabel('Cell Type (Alphabetical)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(base_dir / "cluster_confidence.png", dpi=200, bbox_inches='tight')

def plot_cluster_confidence(summary_df, base_dir: Path):
    clusters = sorted(summary_df['assigned_cluster'].unique())
    fig, ax = plt.subplots(figsize=(max(20, len(clusters) * 0.25), 8))
    
    # TOP PANEL: Max Probability Box Plot
    sns.boxplot(
        data=summary_df, x='assigned_cluster', y='max_probability', 
        order=clusters, ax=ax, color='skyblue', 
        linewidth=1, fliersize=2  # fliersize=2 keeps the min/max dots small so they don't clutter!
    )
    sns.stripplot(
        data=summary_df, x='assigned_cluster', y='max_probability',
        order=clusters, ax=ax, color='black', alpha=0.3, size=2, jitter=True
    )
    ax.set_ylabel('Probability Score', fontsize=14, fontweight='bold')
    ax.set_title(f'Assignment Confidence Inspection ({len(clusters)} Clusters) - Distribution of Cluster Assigned Cells/Nuclei', fontsize=18)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(axis='y', linestyle=':', alpha=0.5)
    ax.set_xlabel('') # Hide the top x-axis label

    # X-Axis formatting
    plt.xticks(rotation=90, fontsize=9)
    plt.xlabel('Cell Type (Alphabetical)', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(base_dir / "cluster_confidence.png", dpi=200, bbox_inches='tight')



def plot_threshold_distribution(summary_df, base_dir: Path):
    thresholds, colors = [0.2, 0.5, 0.7], ['#c7e9b4', '#41b6c4', '#225ea8'] 
    all_clusters = sorted(summary_df['assigned_cluster'].unique())
    fig, axes = plt.subplots(len(thresholds), 1, figsize=(max(20, len(all_clusters) * 0.3), 6 * len(thresholds)), sharex=True)
    
    for i, t in enumerate(thresholds):
        counts = summary_df[summary_df['max_probability'] >= t].groupby('assigned_cluster').size().reindex(all_clusters, fill_value=0)
        axes[i].bar(all_clusters, counts, color=colors[i], alpha=0.8, width=0.8)
        axes[i].set_title(f'Cells with Assignment Probability >= {t}', fontsize=20, fontweight='bold')
        axes[i].set_ylabel('Cell Count', fontsize=14)
        axes[i].grid(axis='y', linestyle='--', alpha=0.3)
        for j, v in enumerate(counts):
            if v > 0: axes[i].text(j, v + (max(counts)*0.01), str(int(v)), rotation=90, ha='center', va='bottom', fontsize=8, alpha=0.7)

    plt.xticks(rotation=90, fontsize=9)
    plt.xlabel('Cell Type (Alphabetical)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(base_dir / "cluster_assignment_distribution.png", dpi=200, bbox_inches='tight')

# --- CLI COMMAND ---
@click.command()
@click.argument("base_dir", type=click.Path(exists=True))
def main(base_dir):
    """Run an audit on the entire dataset (Confidence Distributions)."""
    base_dir_path = Path(base_dir)
    cell_summary = generate_confidence_summary(base_dir_path)
    if cell_summary is not None:
        plot_cluster_confidence(cell_summary, base_dir_path)
        plot_threshold_distribution(cell_summary, base_dir_path)

if __name__ == "__main__":
    main()