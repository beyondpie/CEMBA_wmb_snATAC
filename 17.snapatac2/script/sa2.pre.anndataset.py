# Deprecated
# Now we use sa2.merge.rmdlt.py to get a complete bmat snap file.
# then do subset on it directly.
import os
import sys
import traceback
from pathlib import Path
import logging
from typing import Dict, List
import shutil


from numba.core.errors import NumbaDeprecationWarning 
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import numpy as np
import snapatac2 as sa2
import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils # noqa: E402
from leiden import cemba #pyright: ignore # noqa: E401, E402, F401
from leiden import cal_silhouette #pyright: ignore # noqa: E402, F401
from leiden import umap #pyright: ignore # noqa: E402

logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "cemba.all.anndataset"
)
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

logger.info(f"Get CEMBA all anndata files into AnnDataSet.")


sample2fragment_file = snakemake.input[0]
cemba_all_file = snakemake.output[0]

with open(sample2fragment_file, 'r') as f:
   files = [l.strip() for l in f.readlines()]
samples = [os.path.basename(a).split(".")[0] for a in files]
sample2files = [(s, f) for s, f in zip(samples, files)]
logger.info(f"Load {len(samples)} samples into AnnDataSet.")
sds = sa2.AnnDataSet(
    adatas = sample2files,
    filename = cemba_all_file
)
logger.info("Update obs_names: [sample].[barcode] .")

obs_names: List[str] = [f"{i}.{j}"
              for i, j in zip(sds.obs['sample'].to_list(), sds.obs_names)]
sds.obs_names = obs_names

sds.close()
logger.info(f"CEMBA.all.AnnDataset is saved at: {cemba_all_file}.")

