import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy
import numpy as np
import pandas as pd
import importlib
import snapatac2 as sa2
import gc

# importlib.reload(sa2)
import snapatac2._snapatac2 as internal

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
# importlib.reload(utils) # type: ignore # noqa: F821
from utils import printfn, head # type: ignore # noqa:E402

# * configs
black_list = Path("../meta/mm10.blacklist.bed")
# samples = ['CEMBA171206_3C', 'CEMBA171207_3C',
#            'CEMBA171213_4B', 'CEMBA190423_10A']
samples = ['CEMBA171207_3C',
           'CEMBA171213_4B', 'CEMBA190423_10A']
bmat_dir = "snapatac2_pp_out/bmat"
cl_dir = "snapatac2_pp_out/cl"
if not os.path.exists(cl_dir):
    os.makedirs(cl_dir, exist_ok=True)


# * load multiple bmat and merge them
## read one anndata
sd = sa2.read(Path(f"{bmat_dir}/{samples[0]}.sa2.bmat.hdf5"))
## set AnnDataSet
files = [(s, Path(f"{bmat_dir}/{s}.sa2.bmat.hdf5")) for s in samples]
l1cl_pp_file = f"{cl_dir}/l1.cl.pp.hdf5"
if os.path.exists(l1cl_pp_file):
    os.remove(l1cl_pp_file)
sds = sa2.AnnDataSet(
    adatas = files,
    filename = f"{cl_dir}/l1.cl.pp.hdf5"
)
sds.obs_names = np.array(sds.obs['sample']) + "." + np.array(sds.obs_names)

# * feature filtering
sa2.pp.select_features(
    adata = sds,
    n_features = sds.n_vars,
    filter_lower_quantile = 0.0,
    filter_upper_quantile = 0.005,
    whitelist = None,
    blacklist = black_list,
    max_iter = 1,
    inplace = True
)

# * low-embedding
sa2.tl.spectral(
    adata = sds,
    n_comps = 50,
    features = 'selected',
    random_state = 0,
    sample_size = None,
    sample_method = 'random',
    distance_metric = 'cosine',
    weighted_by_sd = True,
    inplace = True
)

# * UMAP
sa2.tl.umap(
    adata = sds,
    n_comps = 2,
    use_rep = 'X_spectral',
    key_added = 'umap',
    inplace = True
)

# * clustering
sa2.pp.knn(
    adata = sds,
    n_neighbors = 50,
    use_dims = None,
    use_rep = 'X_spectral',
    ## use exact to compare the results
    method = 'hora',
    inplace = True,
    random_state = 0
)
# * resolution selection
sa2.tl.leiden(
    adata = sds,
    resolution = 1.0,
    # same as we used in SnapATAC
    objective_function = 'RBConfiguration',
    min_cluster_size = 5,
    n_iterations = -1,
    random_state = 0,
    key_added = 'leiden',
    use_leidenalg = False,
    inplace = True,
    weighted = False
)

# * UMAP and tsne




