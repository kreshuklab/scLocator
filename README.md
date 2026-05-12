# scLocator
Mapping scRNAseq clusters to images via spatial transcriptomics

## Getting Started

Follow the steps below to install and run scLocator.

### Download the Repository
Download the `scLocator` repository to your local machine. To do this, open terminal and run:
```bash
git clone [https://github.com/kreshuklab/scLocator.git](https://github.com/kreshuklab/scLocator.git)
cd scLocator
```

### Install Dependencies
Create and activate the conda environment used for running the pipeline using:

```
conda env create -f environment.yml
conda activate sclocator
```

## Input Data Requirements
* **scRNA-seq Atlas:** Must be a valid `.h5ad` file with raw or normalized counts in a designated layer, and a metadata column containing cluster annotations.
* **Spatial Data:** Must be a tabular format (`.tsv` or `.csv`) where rows are cells/nuclei and columns are genes. **Crucially, the file must contain a `label_id` column** serving as the unique identifier for each cell/nucleus.

### Run the Pipeline
Once your environment is active and the configuration file is prepared (see details below), the pipeline can be launched from terminal:
```bash
sclocator --config data/experiments/config.yaml
```

## Outputs
After a successful run, the pipeline will generate the following structure in the specified `output_dir`:
* `spatial_expr_df_binarized.tsv` - The joined and binarized spatial expression matrix.
* `cluster_prob/`
  * `cluster_probability.tsv` - The final NxM probability matrix (Cells/Nuclei x Clusters).
  * `probability_heatmap.png` - A visual summary of the assignments.
* *`failed_cell_details.tsv` - (If diagnostics enabled) Cells that failed assignment.*
  * *`generic_cluster_diagnostics.tsv` - (If diagnostics enabled) Cells assigned to noise.*
  * *`gene_likelihood_vectors/` - (If diagnostics enabled) Detailed Bayesian math samples.*
* `expression_profiles/` - Contains the calculated scRNA-seq profiles (`.pkl`) and QC plots.
* `targeted_inspections/` - Created when using the `sclocator-inspect-cell` command.



### Diagnostic Tools
If you want to inspect the results of your mapping, `scLocator` comes with two built-in command-line tools to visualize the assignment math and catch over-clustering or mapping failures.

#### 1. Global Cluster Inspection
This command generates a global summary of your mapping quality. It produces a detailed Box Plot (overlaid with a single-cell scatter plot) showing the distribution of the Best-Match Probability Scores for every single cell type in your atlas.

**Why use this?** To instantly identify clusters that are struggling to map (low probability scores) or to spot hidden bi-modal distributions indicating a cluster might need to be split.

**Usage:**
```bash
sclocator-inspect-clusters <path_to_output_dir>
```
Output: Saves cluster_confidence.png to the target directory.

#### 2. Single-Cell Inspection
If a specific cell is mapped poorly (or simply to see the exact Bayesian math behind an assignment), individual an cell can be inspected by its ID. This generates a detailed heatmap showing exactly how much each gene contributed to (or penalized) the assignment score for the top 10 cluster choices.

Why use this? To figure out exactly which spatial genes are driving the algorithm's decisions, and catch potential errors in spatial gene registration.

**Usage:**
``` bash
sclocator-inspect-cell <nucleus_id> <base_dir> [OPTIONS]
```
**Required Arguments:**
* `nucleus_id` *(Integer/Float)*: The ID of the cell/nucleus you want to inspect.
* `base_dir` *(Path)*: The path to the output directory from your main pipeline run.

**Optional Flags:**
* `top_n` <INT>: The number of top probable clusters to include in the heatmap comparison. (Default: 10)
* `compare_to` <TEXT>: Explicitly force a specific cluster into the comparison, even if it didn't score in the top N. You can use this flag multiple times.
* `gene_id_col` <TEXT>: The column name in your metadata file for raw gene IDs. (Default: "XLOC")
* `gene_name_col` <TEXT>: The column name in your metadata file for readable gene names. (Default: "Platybrowser name")

Output: Saves a targeted heatmap (e.g., cell_8449_heatmap.png) and the raw probability math (cell_8449_vs_targets.tsv) inside the targeted_inspections folder of your output directory.


## Algorithm

The algorithm used in the pipeline to assign cluster probabilities to spatial cells/nuclei is detailed below:

**Given** 
- scRNAseq atlas consisting of $G$ genes and $N_{cells}$ cells split into $N$ clusters
- For gene $k$ number of cells in which this gene is expressed in cluster $i$ is $nonzero_{k,i} \in [0; 1]$

**Input**
- Cell in EM with expression vector $\vec{e}$, for example, $(0, 0, 1, 1, 0)$, so binary measurements "expressed/ not expressed" for a subset of genes of the original atlas. $\vec{e}_g$ is the binary number corresponding to the gene $g$.

**Output**

Probability that given the expression measurements the cell belongs to cluster $i$: $$P(cell \in cluster_i | \vec{e}), i \in \{1 .. N\}$$

**Method**

Using Bayes theorem  
$$P(cell \in cluster_i | \vec{e}) = \frac{P(\vec{e}|cell \in cluster_i)\cdot P(cell \in cluster_i)}{P(\vec{e})}$$

All parts can be easily estimated from the scRNAseq atlas:

$$P(\vec{e}|cell \in cluster_i) = \prod\limits_{g=0}^{G}[\vec{e}_g \cdot nonzero_{g,i} + (1 -\vec{e}_g) \cdot (1 - nonzero_{g,i})]$$

$$P(cell \in cluster_i) = \frac{1}{N_{clusters}}$$
$$P(\vec{e}) = \sum\limits_{k=0}^{N_{clusters}}P(\vec{e}|cell \in cluster_k)\cdot P(cell \in cluster_k)$$


Result is a table of size $N_{cells} \times N_{clusters}$ where each element is a probability that given cell belongs to the given cluster.


To accomodate for the missing signal and autofluorescent areas, two artificial clusters are added:
- **zero**: expression profile with value $0.01$ for each gene, describes cells in which no expression was measured
- **autofluorescence$**:  expression profile with value $0.99$ for each gene, describes cells in which expression appears to be measured for all genes, indicating areas with some artifacts



## Configuration File Parameters

Two example configuration files are provided in the `/tests/` directory: one for high-resolution clustering (small groups) and one for broad cell-type families (large groups).

### General Run Settings (`[run]`)

| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `output_dir` | Path | Directory where the pipeline outputs will be saved. | `"tests/test_output"` |
| `diagnostics` | Boolean | Toggles internal debugging. If `true`, saves data for cells/nuclei that failed the assignment process (e.g., `failed_cell_details.tsv`). | `true` |
| `visualize` | Boolean | Toggles the generation of Quality Control (QC) plots during data preparation. | `true` |
| `vis_name_metadata_column` | String | The metadata column used to display human-readable gene names in plots and outputs. | `"Platybrowser name"` |
| `gene_metadata_path` | Path | Path to the gene metadata file. Used to cross-reference and map genes between the scRNA-seq atlas and spatial datasets. | `"data/scrnaseq_atlas/genes_metadata.tsv"` |
| `use_spatial_weights` | Boolean | Enables spatial specificity weighting. If `true`, scales a gene's likelihood importance based on its spatial distribution, penalizing genes that are ubiquitously expressed. | `true` |

### Reference Atlas Settings (`[atlas]`)

| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `adata_path` | Path | Path to the scRNA-seq reference atlas (`.h5ad` file). | `"data/scrnaseq_atlas/atlas_6dpf.h5ad"` |
| `cluster_column` | String | The metadata column in the AnnData object containing the target cluster annotations. | `"curated_high_res_annotations"` |
| `normalize_expression` | Boolean | Applies standard normalization to the scRNA-seq counts prior to processing. *(Recommended)* | `true` |
| `layer` | String | The name of the data layer in the scRNA-seq atlas to use for expression values. | `"RNA_counts"` |
| `metadata_column` | String | The primary ID column used to join the atlas, spatial data, and metadata lookup file. | `"XLOC"` |
| `bin_thresh` | String/Float | The threshold at which scRNA-seq cluster gene expression is binarized. Accepts `"median"` or a fixed float `(0.0 - 1.0)`. | `"median"` |


### Spatial Datasets (`[spatial_data]`)

This section defines the spatial gene expression modalities used in the pipeline. Because spatial imaging data may come from multiple staining (e.g., HCR-FISH, CISH, IF) and microscopy (e.g., Spinnning Disk, Brightfield, Confocal, Light-sheet) modalities, this section accepts a **list** of datasets. 

You can add as many modalities as you need by starting a new block with a hyphen (`-`). 

**Example Structure:**
```yaml
spatial_data:
  - name: "smfish"
    df_path: "data/smfish/expression_matrix_spots.tsv"
    # ... other parameters ...
  - name: "prospr"
    df_path: "data/smfish/prospr_genes.tsv"
    # ... other parameters ...
```
**Parameters per Dataset:**

| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `name` | String | A unique identifier name for the spatial modality. Used for output folder naming and internal tracking. | `"smfish"` |
| `df_path` | Path | Path to the registered gene expression matrix (TSV/CSV) for this specific modality. | `"data/smfish/matrix.tsv"` |
| `metadata_column` | String | The column name in the master metadata file that matches the gene formatting in this specific expression matrix. | `"Platybrowser name"` |
| `bin_thresh` | Float/Int | The threshold value used to binarize gene expression for this specific modality. Value type depends if data is normalized or raw (if `apply_norm` = true or false). | `20` or `0.05` |
| `apply_norm` | Boolean | Toggles whether the pipeline uses the raw counts (`false`) or the percentile-normalized values (`true`) when applying the binarization threshold. | `false` |
| `pmin/pmax` | Integer | The minimum/maximum percentile bound used if `apply_norm` is `true`. (0-100) | `0/99` |


### Expression Profiles (`[expression_profiles]`)

This section defines the mathematical parameters used when building the cluster expression profiles (the reference probabilities) from the scRNA-seq atlas. These settings help stabilize the Bayesian calculations and handle dropout/zero-count genes.

**Parameters:**

| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `norm_per_gene` | Boolean | Toggles whether to normalize the expression probabilities independently for each gene across all clusters (rather than globally or per-cluster). | `false` |
| `base_prob` | Float | The baseline smoothing probability. Acts as both a floor (for genes with zero expression) and a ceiling (for genes with 100% expression) to prevent absolute 0.0 or 1.0 probabilities, which would break the log-likelihood math. | `0.001` |
| `delta` | Float | A tiny probability assigned to the synthetic "zero" cluster (which catches empty/junk cells). Prevents `log(0)` errors for completely empty expression vectors. | `0.00001` |




