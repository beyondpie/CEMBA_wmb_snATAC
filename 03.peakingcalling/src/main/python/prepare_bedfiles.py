import os
import sys
import logging
from typing import Dict, List
import pandas as pd
import numpy as np
import snapatac2 as sa2

import pyprojroot
proj_dir = str(pyprojroot.here())
pack_dir = f"{proj_dir}/package/python"
sys.path.insert(0, pack_dir)
from mylog import StreamToLogger, set_file_logger
from utils import complete_values

# * snakemake
all_snap2_fnm = snakemake.input["snap2"]
b2c_fnm = snakemake.input["b2c"]
clusters = snakemake.params["clusters"]
log_fnm = snakemake.log[0]

outdir = snakemake.params["outdir"]
prefix = snakemake.params["prefix"]
suffix = snakemake.params["suffix"]
bedtag = snakemake.wildcards["bedtag"]

myseed = snakemake.params["seed"]
debug = snakemake.params["debug"]


# * debug
if debug > 0:
    print("Under debug mode...")
    # rsc_dir = os.path.join(proj_dir, "18.snap2_peakcalling",
    #                     "src/main/resource")
    # all_snap2_fnm = os.path.join(rsc_dir, "cemba.snap2.with.fragment.hdf5")
    # b2c_fnm = os.path.join(rsc_dir, "neuron_barcode2L4p2rep_cca-k50-cl_v1.csv")
    # pL4_fnm = os.path.join(rsc_dir, "neuron_L4pc2sizes_cca-k50-cl_v1.csv")
    # pL4meta = pd.read_csv(pL4_fnm, sep = ',', header = None,
    #                     names = ["cluster", "size", "early_size", "later_size"],
    #                     index_col = None)
    # clusters: List[str] = pL4meta['cluster'][[2,5]]
    clusters: List[str] = ['2-25-2-1']
    # outdir = os.path.join(proj_dir, "18.snap2_peakcalling",
    #                     "out")
    # log_fnm = os.path.join(outdir, "log", "tscc_test_bedfile.log")
    # prefix = "neuron."
    # suffix = ".bed.gz"
    # myseed = 0
else:
    print("Under normal mode ...")

# * log
logger = set_file_logger(log_fnm, name = "prepare_bedfiles")
sys.stdout = StreamToLogger(logger=logger, level=logging.INFO)
sys.stderr = StreamToLogger(logger=logger, level=logging.ERROR)

# * main
if isinstance(clusters, str):
    logger.info(f"clusters is one str: {clusters}")
    logger.info("make it a list.")
    clusters = [clusters]
logger.info(f"# of clusters: {len(clusters)}")

logger.info(f"load snap2 object from: {all_snap2_fnm}.")
## read only without modificaiton
all_snap2 = sa2.read_dataset(all_snap2_fnm, mode = 'r')
all_barcodes: List[str] = all_snap2.obs_names

logger.info("get groupby column for all_snap2.")
logger.info(f"load barcode2cluster from: {b2c_fnm}.")
b2c: pd.DataFrame = pd.read_csv(b2c_fnm,
                                sep = ',',
                                header = None,
                                names = ["barcode", "cluster"],
                                index_col = None)


logger.info("get all fragments for clusters.")

b2c_dict = dict(zip(b2c.barcode, b2c.cluster))

all_clusters: List[str] = complete_values(key2value = b2c_dict,
                               all_keys = all_barcodes,
                               default_value = "NA")

if (len(clusters) < 2) and ( clusters[0] == "NA"):
    logger.info("clusters is NA, will use all the ones in b2c file.")
    clusters: List[str] = b2c.cluster.unique().tolist()
    logger.info(f"Now # of clusters: {len(clusters)}.")
    
    
sa2.ex.export_bed(adata = all_snap2,
                  groupby = all_clusters,
                  selections = clusters,
                  out_dir = outdir,
                  prefix = prefix,
                  suffix = suffix)
logger.info("Done.")

