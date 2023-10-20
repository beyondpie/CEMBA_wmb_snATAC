from os.path import exists
import pandas as pd

# 1.23.0 since
# - numba(cell oracle needs) <=1.24, but
# - seaborn != 1.24.0

import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
import pyarrow.parquet as pq
import pyprojroot
import matplotlib.pyplot as plt

proj_root = pyprojroot.here()
print(proj_root)
co_base_sc_dir = os.path.join(proj_root, "22.sa2GRN",
                              "out", "GRN",
                              "baseGRN_subclass")

out_dir = os.path.join(proj_root, "22.sa2GRN",
                       "out", "powerlaw")
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok = True)

sc = "Astro-TE_NN"
link_sc = co.utility.load_hdf5(
    os.path.join(co_base_sc_dir, f"GRN.{sc}.celloracle.links"))

link_sc
link_sc.cluster

# use default parameter for filtering
# p=0.001, weight="coef_abs", threshold_number=10000
link_sc.filter_links()
link_sc.plot_degree_distributions(plot_model=True,
                                  save = os.path.join(out_dir, sc))
