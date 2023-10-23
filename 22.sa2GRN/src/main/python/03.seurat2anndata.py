import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
from anndata import AnnData
from scipy.io import mmread
from scipy.sparse import csr_matrix
import pyarrow.parquet as pq
import pyprojroot
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

proj_dir = pyprojroot.here()
packdir = f"{proj_dir}/package/python"
sys.path.insert(0, packdir)
import importlib
# importlib.reload(utils)
import utils
# importlib.reload(mycelloracle)
import mycelloracle

# * load mtx to anndata
# Ref: https://github.com/morris-lab/CellOracle/blob/master/celloracle/data_conversion/process_seurat_object.py
# Ref: https://github.com/scverse/anndata/blob/main/anndata/_io/read.py#L303-L321

mtx_dir = os.path.join(proj_dir, "22.sa2GRN", "src/main/resource",
                       "sa2.allen.logCPM.vf3281.ds1000")
X = mmread(os.path.join(
    mtx_dir,"sa2.allen.logCPM.vf3281.ds1000.gbyc.mtx")).astype("float32")
X = csr_matrix(X)
X_cbyg = csr_matrix.transpose(X)
mm = AnnData(X_cbyg)
cell_ids = pd.read_csv(os.path.join(mtx_dir, "barcode.txt"),
                       header = None).values[:,0]
gene_ids = pd.read_csv(os.path.join(mtx_dir, "gene.txt"),
                       header = None).values[:, 0]
meta_info = pd.read_csv(os.path.join(mtx_dir, "meta.data.csv"))

mm = AnnData(X = X_cbyg,
             obs = meta_info.reindex(cell_ids),
             var = pd.DataFrame(index = gene_ids))
# ... storing 'orig.ident' as categorical
# ... storing 'modality' as categorical
# ... storing 'library_prep' as categorical
# ... storing 'roi' as categorical
# ... storing 'method' as categorical
# ... storing 'sex' as categorical
# ... storing 'age' as categorical
# ... storing 'medical_conditions' as categorical
# ... storing 'subclass' as categorical
mm.write_h5ad(
    filename = os.path.join(mtx_dir,
        "sa2.allen.logCPM.vf3281.ds1000.h5ad"))

# * save each subclass for later calculation
out_dir = os.path.join(proj_dir, "22.sa2GRN", "src/main/resource",
                       "sa2.allen.logCPM.vf3281.ds1000.subclass.specific")
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok = True)
subclasses = mm.obs['subclass'].unique()
for sc in subclasses:
    print(f"Get subclass {sc}")
    adata_local = mm[mm.obs['subclass'].isin([sc])]
    adata_local.write(os.path.join(out_dir, f"sa2.allen.subclass.{sc}.ann.hdf5"))

# * select subclass with at least some cells
# min_ncell = 50
# nstat_subclass = mm.obs.groupby(['subclass']).size()

# cl2 = concatl2_to_n[concatl2_to_n >= 50].index

# L3IntMeta = pd.read_csv("../meta/IntL3SumMeta.csv", sep = ",",
#                         header = 0, index_col = None)
# mappedl2 = L3IntMeta['L3AllenAnnot'].apply(utils.concat_allen_name).unique()

# myl2 = pd.DataFrame(cl2[cl2.isin(mappedl2)].to_list())
# myl2.to_csv("allen.subclass.chosen.txt", sep = ",", index = False, header = False)

