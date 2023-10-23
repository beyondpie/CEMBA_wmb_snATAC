import os
import sys
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import snapatac2 as sa2

import pyprojroot

proj_dir = str(pyprojroot.here())
pack_dir = f"{proj_dir}/package/python"
sys.path.insert(0, pack_dir)
from mylog import StreamToLogger, set_file_logger
from mysnapatac2 import modify_obs_name

# * snakemake
log_fnm: str =snakemake.log[0]
bedfnm: str = snakemake.input["bedfnm"]
snap2_dir: str = snakemake.params["snap2_dir"]
suffix: str = snakemake.params["suffix"]
out_dir: str = snakemake.params["out_dir"]
out_suffix: str = snakemake.params["out_suffix"]
sample: str = snakemake.wildcards["s"]

# * debug
# log_fnm = "test_pmat.log"
# work_dir = os.path.join(proj_dir, "18.snap2_peakcalling")
# bedfnm = os.path.join(proj_dir, "18.snap2_peakcalling",
#                       "out/tscc", "mba.whole.union.peak.srt.bed")
# with open(f"{proj_dir}/17.snapatac2/meta/mba.whole.sample.lst", 'r') as f:
#     samples = [l.strip() for l in f.readlines()]

# snap2_dir = os.path.join(proj_dir,
#                          "17.snapatac2",
#                          "sa2_qc_dlt/rm_dlt")
# out_dir = "."
# suffix = "_rm_dlt.h5ad"
# out_suffix = "_union_pmat.h5ad"
# sample = samples[0]

# * log
logger = set_file_logger(log_fnm, name="sa2.get_full_snap2")
sys.stdout = StreamToLogger(logger=logger, level=logging.INFO)
sys.stderr = StreamToLogger(logger=logger, level=logging.ERROR)

# * main
# ** load snap2
logger.info(f"Loading snap2 for {sample}...")
snap_file = os.path.join(snap2_dir, f"{sample}{suffix}")
if not os.path.exists(snap_file):
    raise FileNotFoundError(f"{snap_file} not found.")
snap2: sa2.AnnData = sa2.read(snap_file, backed = 'r')

# ** read bed file
logger.info("Reading bed file...")
if not os.path.exists(bedfnm):
    raise FileNotFoundError(f"{bedfnm} not found.")


# ** get pmat
logger.info(f"Getting pmat for {sample}...")
outfnm = os.path.join(out_dir, f"{sample}{out_suffix}")
sa2.pp.make_peak_matrix(adata = snap2,
                        inplace = False,
                        file = outfnm,
                        backend = "hdf5",
                        peak_file = bedfnm,
                        chunk_size = 10000,
                        use_x = False)
snap2.close()
logger.info(f"Done.")


