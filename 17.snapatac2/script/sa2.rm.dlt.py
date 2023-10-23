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
    name = "sa2.embed"
)
# logger = utils.set_file_logger( #pyright: ignore
#     fnm = "test_qc_dlt.log", #pyright: ignore # noqa: F821
#     name = "sa2.embed"
# )
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

snap_file = snakemake.input["qc_dlt_file"][0] #pyright: ignore # noqa: F821
out_file = snakemake.output["snap_file"][0] #pyright: ignore # noqa: F821

logger.info(f"Load snap file {snap_file}")
snap = sa2.read(snap_file, backed = 'r')

barcodes: np.ndarray = np.array([snap.obs_names])
dlt_probs = snap.obs['doublet_probability'].to_numpy()
slt_index = (dlt_probs <= 0.5).tolist()

r = snap.subset(obs_indices = slt_index, out = out_file)
r.close()
snap.close()

