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

if __name__ == "__main__":
    main()
