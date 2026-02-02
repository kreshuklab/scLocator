import logging


def validate_spatial_source_config(spatial_data_source):
    logging.info("Validate config for spatial data")
    try:
        if not type(spatial_data_source["name"]) is str:
            raise ValueError

    except KeyError:
        logging.error("name is obligatory part of the spatial data source config")
        raise

    except ValueError:
        logging.error(f"name is {spatial_data_source['name']}, should be a string")
        raise

    try:
        if not type(spatial_data_source["df_path"]) is str:
            raise ValueError

    except KeyError:
        logging.error("df_path is obligatory part of the spatial data source config")
        raise

    except ValueError:
        logging.error(
            f"df_path is {spatial_data_source['df_path']}, should be a string"
        )
        raise

    try:
        bin_thresh = float(spatial_data_source["bin_thresh"])

    except KeyError:
        logging.error("bin_thresh is obligatory part of the spatial data source config")
        raise

    except ValueError:
        logging.error(
            f"bin_thresh is {spatial_data_source['bin_thresh']}, should be a number"
        )
        raise

    spatial_data_source["bin_thresh"] = bin_thresh

    try:
        apply_norm = spatial_data_source["apply__norm"]

    except KeyError:
        spatial_data_source["apply_norm"] = False
        spatial_data_source["pmin"] = 0
        spatial_data_source["pmax"] = 100

    try:
        pmin = float(spatial_data_source["pmin"])
        if pmin < 0 or pmin > 100:
            raise ValueError

    except KeyError:
        pmin = 0

    except ValueError:
        logging.error(
            f"pmin is {spatial_data_source['pmin']}, should be a number between 0 and 100"
        )
        raise

    spatial_data_source["pmin"] = pmin

    try:
        pmax = float(spatial_data_source["pmin"])
        if pmax < 0 or pmax > 100:
            raise ValueError

    except KeyError:
        pmax = 0

    except ValueError:
        logging.error(
            f"pmax is {spatial_data_source['pmax']}, should be a number between 0 and 100"
        )
        raise

    spatial_data_source["pmin"] = pmax

    return spatial_data_source


def validate_atlas_config(atlas_config):
    logging.info("Validate config for scRNAseq atlas")
    try:
        if not type(atlas_config["adata_path"]) is str:
            raise ValueError

    except KeyError:
        logging.error("adata_path is an obligatory parameter")
        raise

    except ValueError:
        logging.error(f"adata_path is {atlas_config['adata_path']}, should be a string")
        raise

    try:
        if not type(atlas_config["cluster_column"]) is str:
            raise ValueError

    except KeyError:
        logging.error("cluster_column is an obligatory parameter")
        raise

    except ValueError:
        logging.error(
            f"cluster_column is {atlas_config['cluster_column']}, should be a string"
        )
        raise

    try:
        if not type(atlas_config["metadata_column"]) is str:
            raise ValueError

    except KeyError:
        logging.error("metadata_column is an obligatory parameter")
        raise

    except ValueError:
        logging.error(
            f"metadata_column is {atlas_config['metadata_column']}, should be a string"
        )
        raise

    try:
        if not type(atlas_config["metadata_column"]) is str:
            raise ValueError

    except KeyError:
        logging.error("metadata_column is an obligatory parameter")
        raise

    except ValueError:
        logging.error(
            f"metadata_column is {atlas_config['metadata_column']}, should be a string"
        )
        raise

    try:
        normalize_expression = atlas_config["normalize_expression"]

    except KeyError:
        atlas_config["normalize_expression"] = False

    try:
        layer = atlas_config["layer"]

    except KeyError:
        atlas_config["layer"] = "X"

    try:
        layer = atlas_config["bin_thresh"]

    except KeyError:
        atlas_config["bin_thresh"] = "median"

    return atlas_config
