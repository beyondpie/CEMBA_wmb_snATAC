import os
import sys
import importlib
import pandas as pd
import anndata as ad

from pyprojroot import here
proj_root: str = str(here())
sys.path.insert(0, os.path.join(proj_root, "package/python"))
from utils import set_file_logger
import cembav2env
importlib.reload(cembav2env)

# * configs
rsc_dir = os.path.join(proj_root,
                       "19.snap2_integration/"
                       "src/main/resource")
allen = cembav2env.Allen()
sa2atac = cembav2env.Sa2ATAC()

# * load reduced ann data
reduced_allen_dir = os.path.join(rsc_dir, "norawdata_allen")
neuron_allen_fnm = os.path.join(reduced_allen_dir,
                            "neuron_male_allen_ann_noraw.h5ad")
neuron_allen = ad.read_h5ad(
    filename=neuron_allen_fnm, backed="r")

nn_allen_fnm = os.path.join(reduced_allen_dir,
                            "nn_male_allen_ann_noraw.h5ad")
nn_allen = ad.read_h5ad(
    filename=nn_allen_fnm, backed="r")

