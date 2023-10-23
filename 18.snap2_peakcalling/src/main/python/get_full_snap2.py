"""
Get snap2 file with all the samples.

"""
import os
import sys
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import snapatac2 as sa2

import pyprojroot

proj_dir = str(pyprojroot.here())
pack_dir = f"{proj_dir}/package/python"
sys.path.insert(0, pack_dir)
from mylog import StreamToLogger, set_file_logger
from mysnapatac2 import modify_obs_name

# * snakemake
# log_fnm=snakemake.log[0]
# snap2_files: List[str] = snakemake.input["snap2_files"]
# out_snap2: str = snakemake.output["snap2"]

# * log
log_dir = os.path.join(proj_dir, "18.snap2_peakcalling", "out/log")
os.makedirs(log_dir, exist_ok=True)
log_fnm = os.path.join(log_dir, "get_full_snap2.log")
logger = set_file_logger(log_fnm, name="sa2.get_full_snap2")
sys.stdout = StreamToLogger(logger=logger, level=logging.INFO)
sys.stderr = StreamToLogger(logger=logger, level=logging.ERROR)

# * meta
rdir = os.path.join(proj_dir, "18.snap2_peakcalling", "src/main/resource")
with open(f"{rdir}/mba.whole.sample.lst", "r") as f:
    samples = [l.strip() for l in f.readlines()]
out_snap2_fnm = f"{rdir}/cemba.snap2.with.fragment.hdf5"

snap2s_dir = os.path.join(proj_dir, "17.snapatac2", "sa2_qc_dlt", "rm_dlt")
snap2_files = [f"{snap2s_dir}/{s}_rm_dlt.h5ad" for s in samples]

logger.info(f"In total, {len(snap2_files)} are inputed.")

fnms = [os.path.basename(v) for v in snap2_files]

logger.info(f"get full snap anndataset object to: {out_snap2_fnm}")
sample2files = [(s, f) for s, f in zip(samples, snap2_files)]
sds = sa2.AnnDataSet(adatas=sample2files,
                        filename=out_snap2_fnm, add_key="sample")
new_obs_names = modify_obs_name(sds, obs_key="sample")
sds.obs_names = new_obs_names
sds.close()
logger.info("Done.")
