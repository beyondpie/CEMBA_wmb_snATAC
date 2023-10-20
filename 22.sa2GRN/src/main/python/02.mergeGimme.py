import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
import pyarrow.parquet as pq
import pyprojroot
packdir = f"{pyprojroot.here()}/package/python"
sys.path.insert(0, packdir)
import utils
import mycelloracle

# * function
def load_tfidf(f) -> pd.DataFrame:
    import pyarrow.parquet as pq
    r = pq.read_table(f)
    return r.to_pandas()

# * load groups
proj_dir = str(pyprojroot.here())
group_file = os.path.join(proj_dir, "meta", "sa2.subclass.srt.txt")
with open(group_file) as f:
    lines = f.readlines()
    groups = [l.strip() for l in lines if len(l.strip()) > 1]
    groups = [g for g in groups if (not "Hypendymal_NN" in g)]

with open("pdc.suffix.txt", mode = "r") as f:
    parts = [l.strip() for l in f.readlines() if len(l) > 1]

# * check gimmemotifs
tfscan_dir = os.path.join(proj_dir, "22.sa2GRN", "out/tfscan")
prefix = "sa2subclass"

# about 100G RAM
tfidfs = [load_tfidf(f"{tfscan_dir}/{prefix}.{i}.pdc.baseGRN.df.parquet")
          for i in groups]
# another 100G RAM
r = pd.concat(tfidfs, axis = 0, ignore_index= True)
r = r.fillna(int(0))
del tfidfs
# 520M
r.to_parquet(f"{tfscan_dir}/{prefix}.all.baseGRN.df.parquet")

# stat
tfCount = r.sum(axis = 0, numeric_only = True)
peakCount = r.sum(axis = 1, numeric_only = True)





