## NOTE: this file is quite similar with sa2.merge.rmdlt.py
##  only diff: get the sample names
## TODO: merge sa2.merge.gmat and sa2.merge.rmdlt.py

import os
import sys
import traceback
from pathlib import Path
from typing import Dict

import numpy as np
import snapatac2 as sa2

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils #pyright: ignore # noqa: F401, E402

# * log
logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "sa2.merge.gmat"
)
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(
                             exc_type, exc_value, exc_traceback)
                         ]))
# Install exception handler
sys.excepthook = handle_exception

snap_files = snakemake.input["snap_files"]
out_snap = snakemake.output["merge_snap"]
tmp_snap = os.path.join(os.path.dirname(out_snap), "tmp.merge.snap.h5ad")
genome = snakemake.params['genome']

logger.info(f"In total, {len(snap_files)} are detected.")

fnms = [os.path.basename(v) for v in snap_files]
samples = [a.replace(f"_{genome}_gmat.h5ad", "") for a in fnms]

sample2files = [(s, f) for s, f in zip(samples, snap_files)]
logger.info(f"Create AnnDataSet to tmp file: f{tmp_snap}")
sds = sa2.AnnDataSet(
    adatas = sample2files,
    filename = tmp_snap,
    add_key = 'sample'
)

logger.info(f"AnnDataSet to AnnData: {out_snap}")

snap = sds.to_adata(file = out_snap, copy_x = True)
new_obs_names = utils.modify_obs_name(snap, obs_key = "sample")
snap.obs_names = new_obs_names

snap.close()
sds.close()
logger.info(f"Delete tmp file: {tmp_snap}.")
os.remove(tmp_snap)
logger.info("Done")
