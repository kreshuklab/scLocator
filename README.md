# scLocator
Mapping scRNAseq clusters to images via spatial transcriptomics

## Installation
```
conda env create -f environment.yml
```

## Algorithm

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


## TODO: explanation of parameters in the config




Create environment using `environment.yml`

Run:

``` python experiments/run_mapping.py --scrnaseq_atlas data/scrnaseq_atlas/updated_6dpf_atlas.h5ad --gene_metadata data/scrnaseq_atlas/genes_metadata.tsv --cluster_column broad.types --smfish_data data/smfish/expression_matrix_aggregated_confocal_fixed.tsv --prospr_data data/smfish/prospr_genes_nuclei.tsv --smfish_thresh 20 --prospr_thresh 0.1 --visualize --run_diagnostic --profile_constant data/experiments/large_groups/ data/experiments/large_groups/cluster_mapping.log```

Here:

`data/scrnaseq_atlas/updated_6dpf_atlas.h5ad` - scRNAseq atlas object

`data/scrnaseq_atlas/gense_metadata.tsv` - metadata table to convert gene names to XLOC (gene IDs in scRNAseq object)

`broad.types` - column in the scRNAseq atalas `obs` to use as cluster names, also made to run 'curated_high_res_annotations' (small_groups)

`data/smfish/expression_matrix_aggregated_confocal_fixed.tsv` - output of the registration pipeline

`data/experiments/large_groups/` - results (plots and tables) are saved in this directory

`data/experiments/large_groups/cluster_mapping.log` - log with some QC information

Optional arguments:
`--smfish_thresh` - defines the thresholds for binarizing HCR stained smfish data, in units of spots, so a threshold of 20 would binarize all nuclei/cells which express more than 20 spotiflow recognized HCR spots
`--prospr_thresh` - defines the thresholds for binarizing prospr data, which comes from WMISH staining and has already been normalized on expression from 0-1
`--visualize` - if passed plots all visualizations
`--run_diagnostic` - if passed outputs diagnostic tables helpful for figuring out certain nuclei/cell assignment probabilities, as well as smfish & prospr gene integration
`--profile_constant` - defines a constant base value added to scRNAseq profiles for each cluster in order to stop probabilities of 0 or 1 of occurring, default is 0.001
 
