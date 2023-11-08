import os
import sys
import math
import pickle
from typing import Dict, List
from dataclasses import dataclass, field
import numpy as np

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
from leiden import LeidenSum, ScatterPlot
from leiden import draw_umap
from leiden import init_LeidenSum_from_file
from colors import SnapATACPalette
import snapatac2 as sa2

sa2_dir = os.path.join(
    "/projects/ps-renlab/szu/projects/CEMBA2",
    "17.snapatac2"
)

sa2L1_fnm = os.path.join(
    sa2_dir, "L1_encoder",
    "nfeat-all_nsample-all_nc50_0_mult.h5ad")

snapL1 = sa2.read(
    filename = sa2L1_fnm, backed = 'r')

embed_mat: np.ndarray = snapL1.obsm['X_spectral']
barcodes: List[str] = snapL1.obs_names
out_dir = os.path.join(
    "/projects/ps-renlab2/szu/projects/CEMBA2",
    "17.snapatac2", "resource", "sa2L1sum"
)
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok = True)

np.savetxt(
    os.path.join(out_dir, "sa2.L1.embed_mat.csv"),
    embed_mat, delimiter = ',')

# save barcodes to txt
with open(os.path.join(
        out_dir, "sa2.L1.barcodes.txt"), 'w') as f:
    for bc in barcodes:
        f.write(bc + '\n')

# save umap
umap: np.ndarray = snapL1.obsm["X_umap"]
np.savetxt(os.path.join(out_dir, "sa2.L1.umap_ab_spectral.csv"),
           umap, delimiter = ',')

# calculate umap with default parameter
from numba.core.errors import NumbaDeprecationWarning 
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
from leiden import umap  # pyright: ignore # noqa: E402

# it will automatically use <20 CPUs
# and cost <30G RAM
# in encoder
umap_default : np.ndarray = umap(
    adata = snapL1,
    use_rep = "X_spectral",
    inplace = False,
    a = None,
    b = None,
    init='spectral'
)

np.savetxt(os.path.join(out_dir, "sa2.L1.umap_default.csv"),
           umap_default, delimiter = ',')
snapL1.close()
