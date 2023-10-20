import sys
import os
import importlib
import anndata as ad
import pandas as pd
import numpy as np
from pyprojroot import here

proj_root = str(here())
sys.path.insert(0, os.path.join(proj_root, "package/python"))
from utils import set_file_logger
import cembav2env

importlib.reload(cembav2env)

rsc_dir = os.path.join(proj_root, "19.snap2_integration",
                       "src/test/resource")
out_allen_dir = os.path.join(rsc_dir, "noraw_allen")
os.makedirs(out_allen_dir, exist_ok = True)
out_atac_dir = os.path.join(rsc_dir, "noraw_atac")
os.makedirs(out_atac_dir, exist_ok = True)

allen = cembav2env.Allen()
sa2atac = cembav2env.Sa2ATAC()

nn_allen = allen.read_nn_10xv3_lognorm_ann()

rdm_nn_allen = nn_allen[np.random.choice(nn_allen.n_obs, 5000, replace = False)]

rdm_nn_allen.write(
    os.path.join(out_allen_dir, "nn_male_10xv3_ann_noraw.h5ad"))

del nn_allen

nn_atac = sa2atac.read_nn_gmat_lognorm_ann()
rdm_nn_atac = nn_atac[np.random.choice(nn_atac.n_obs, 5000, replace = False)]

rdm_nn_atac.write(
    os.path.join(out_atac_dir, "nn_gmat_atac_ann.h5ad"))
del nn_atac

