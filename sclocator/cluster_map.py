import yaml

# from sclocator.preprocess_atlas import get_atlas
from pathlib import Path
import logging
import sys
import pandas as pd
import numpy as np
import scanpy as sc
from sclocator.preprocess_spatial_data import (
    preprocess_spatial_expr_df,
    join_spatial_expr_dfs,
    harmonize_gene_names,
    calculate_spatial_specificity_weights,
    calculate_gene_specificity_weights
)
from sclocator.config_validation import (
    validate_spatial_source_config,
    validate_atlas_config,
)
from sclocator.preprocess_atlas import get_atlas, get_cluster_profiles, visualize_profiles
from sclocator.bayesian_prob import assign_clusters



class ClusterMapper:
    def __init__(self, config_path):
        self.analysis_setup(config_path)

        # self.atlas = get_atlas(**config["atlas"])

    def analysis_setup(self, config_path):
        with open(config_path, "r") as yamlfile:
            config = yaml.load(yamlfile, Loader=yaml.FullLoader)

        self.config = config
        self.diagnostics = config["run"]["diagnostics"]
        self.visualize = config["run"]["visualize"]
        self.vis_name_metadata_column = config["run"]["vis_name_metadata_column"]
        pd.set_option("display.max_columns", None)
        pd.set_option("display.max_rows", None)

        self.output_dir = Path(config["run"]["output_dir"])
        self.output_dir.mkdir(exist_ok=True)

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler(self.output_dir / "sclocator.log", mode="w"),
                logging.StreamHandler(sys.stdout),
            ],
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        logging.info(f"Finished setup with config at {config_path}")

    def run_analysis(self):
        logging.info(f"Start analysis")
        self.preprocess_data()
        self.create_expression_profiles()
        self.get_cluster_prob()

    def preprocess_data(self):
        self.metadata = pd.read_csv(self.config["run"]["gene_metadata_path"], sep="\t")

        spatial_dfs = []
        df_names = []
        for spatial_data_source in self.config["spatial_data"]:
            spatial_data_source = validate_spatial_source_config(spatial_data_source)

            spatial_data_source["visualize"] = self.visualize
            spatial_data_source["output_dir"] = (
                self.output_dir / spatial_data_source["name"]
            )
            spatial_data_source["output_dir"].mkdir(exist_ok=True)
            logging.info(spatial_data_source["name"])
            expr_df = preprocess_spatial_expr_df(**spatial_data_source)
            expr_df = harmonize_gene_names(
                self.metadata,
                expr_df,
                self.config["atlas"]["metadata_column"],
                spatial_data_source["metadata_column"],
            )
            spatial_dfs.append(expr_df)
            df_names.append(spatial_data_source["name"])

        self.spatial_expr_df = join_spatial_expr_dfs(
            spatial_dfs,
            df_names,
            self.metadata,
            self.config["atlas"]["metadata_column"],
            self.config["run"]["vis_name_metadata_column"],
        )

        self.spatial_expr_df.to_csv(
            self.output_dir / "spatial_expr_df_binarized.tsv",
            sep="\t",
            index=True,
            index_label="label_id",
        )

        self.spatial_gene_subset = self.spatial_expr_df.columns

        self.atlas_config = validate_atlas_config(self.config["atlas"])
        self.atlas_config["visualize"] = self.visualize
        self.atlas_config["output_dir"] = self.output_dir / "atlas"
        self.atlas_config["output_dir"].mkdir(exist_ok=True)
        self.atlas_config["metadata"] = self.metadata
        self.atlas_config["gene_list"] = self.spatial_gene_subset
        self.atlas_config["vis_name_metadata_column"] = self.config["run"][
            "vis_name_metadata_column"
        ]
        self.atlas = get_atlas(**self.atlas_config)
        # self.atlas = sc.read_h5ad(self.atlas_config["output_dir"] / "scrnaseq_bin_expr.h5ad")

        self.spatial_expr_df = self.spatial_expr_df[self.atlas.var.index]

        self.gene_list = self.atlas.var.index

    def create_expression_profiles(self):
        output_dir = self.output_dir / "expression_profiles"
        output_dir.mkdir(exist_ok=True)
        self.profiles = get_cluster_profiles(
            self.atlas,
            self.config["atlas"]["cluster_column"],
            output_dir,
            self.config["expression_profiles"]["norm_per_gene"],
            self.config["expression_profiles"]["base_prob"],
            self.config["expression_profiles"]["delta"],
        )

        visualize_profiles(self.atlas.var.vis_name, self.profiles, output_dir)

    def get_cluster_prob(self):
        output_dir_cluster_prob = self.output_dir / "cluster_prob"
        output_dir_cluster_prob.mkdir(exist_ok=True)

        if self.config["run"]["use_spatial_weights"]:
            spatial_weights = calculate_spatial_specificity_weights(self.spatial_expr_df)
            sc_weights = calculate_gene_specificity_weights(self.profiles, self.gene_list)
            gene_weights = spatial_weights
        else:
            gene_weights = np.ones(len(self.gene_list))
        assign_clusters(self.spatial_expr_df, self.gene_list, self.profiles, output_dir_cluster_prob, self.atlas_config["cluster_column"], gene_weights, False, output_dir_cluster_prob)


    def visualize_result(self):
        pass

    def write_result(self):
        pass
