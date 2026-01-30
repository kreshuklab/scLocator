import os
import sys
import click
import logging
import numpy as np
import matplotlib.pyplot as plt

from sclocator.cluster_map import ClusterMapper



# @click.option("-fi", "--fixed_path", required=True, help="Fixed input .n5 file")
# @click.option("-fk", "--fixed_key", required=True, help="Fixed input key")
# @click.option("-mi", "--moving_path", required=True, help="Moving input .n5 file")
# @click.option("-mk", "--moving_key", required=True, help="Moving input key")
# @click.option("-o", "--output_dir", required=True, help="Output directory")
# @click.option("-ok", "--output_key", required=True, help="Output key (same in both n5)")
# 
@click.command()
@click.option("-c", "--config", required=True, help="Path to the config file", type=str)
def main(config):

    obj = ClusterMapper(config)

    obj.run_analysis()
    # print(config)
    # logging.basicConfig(
    #     level=logging.INFO,
    #     format="%(asctime)s [%(levelname)s] %(message)s",
    #     handlers=[
    #         logging.FileHandler(f"log.log", mode="w"),
    #         logging.StreamHandler(sys.stdout),
    #     ],
    #     datefmt="%Y-%m-%d %H:%M:%S",
    # )
    
    # logging.info("Reading fixed image")
    # fixed_img = read_volume(fixed_path, fixed_key)
    # logging.info(f"Fixed image shape: {fixed_img.shape}, dtype {fixed_img.dtype}")

    # logging.info("Reading moving image")
    # moving_img = read_volume(moving_path, moving_key)
    # logging.info(f"Moving image shape: {moving_img.shape}, dtype {moving_img.dtype}")

    # fixed_prealigned, moving_prealigned, prealignment_transform = run_prealignment(
    #     fixed_img,
    #     moving_img,
    #     output_dir
    # )

    # logging.info("Save prealigned fixed image")
    # fixed_attributes = dict(get_attrs(fixed_path, fixed_key))
    # write_volume(
    #     f=fixed_path,
    #     arr=fixed_prealigned,
    #     key=output_key,
    #     attrs=fixed_attributes
    # )

    # logging.info("Save prealigned moving image")
    # moving_attributes = dict(get_attrs(moving_path, moving_key))
    # write_volume(
    #     f=moving_path,
    #     arr=moving_prealigned,
    #     key=output_key,
    #     attrs=moving_attributes
    # )

    # logging.info("Save prealignment tranform")

    # print(prealignment_transform)

    # write_transform_dict(prealignment_transform, output_transform_path)


if __name__ == "__main__":
    main()
