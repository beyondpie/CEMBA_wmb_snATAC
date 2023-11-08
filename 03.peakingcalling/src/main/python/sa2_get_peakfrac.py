import os
import sys

import logging
from pathlib import Path
from typing import Dict, List, Tuple
from multiprocessing import Pool
import pickle

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import snapatac2 as sa2
import anndata
# import scanpy as sc

import pyprojroot
import pyreadr

proj_dir = str(pyprojroot.here())
pack_dir = f"{proj_dir}/package/python"
sys.path.insert(0, pack_dir)
from mylog import StreamToLogger, set_file_logger
from mysnapatac2 import modify_obs_name
from cembav2env import rename_allensubclass, read_bed4

# * meta
## lastest metadata version: v9.6
atacmetafnm = os.path.join(proj_dir,
                           "supple.02.annotation.all.in.one.meta",
                           "mba.whole.cell.meta.v9.6.rds")
atacmeta: pd.DataFrame = pyreadr.read_r(atacmetafnm)[None]
atacmeta.set_index("barcode2", inplace = True)

sa2pmatd = os.path.join(proj_dir, "18.snap2_peakcalling",
                        "out/tscc/sa2pmat")
samples: List[str] = atacmeta["sample"].unique().tolist()

sa2pmatd = os.path.join(proj_dir, "18.snap2_peakcalling",
                        "out/tscc/sa2pmat")

# * get AnnDataSet and AnnData for union peak
sa2pmatd2 = os.path.join(sa2pmatd, "union_pmat")
suffix = "_union_pmat.h5ad"
sample2files: List[Tuple[str, str]] = [
    (sample, os.path.join(sa2pmatd2, f"{sample}{suffix}")) for sample in samples
]
outsdsfnm = os.path.join(proj_dir, "18.snap2_peakcalling",
                         "out/scfilter/sa2pmat.union_pmat.h5ad")
if not os.path.exists(os.path.dirname(outsdsfnm)):
    os.makedirs(os.path.dirname(outsdsfnm), exist_ok = True)
    
sds: sa2.AnnDataSet = sa2.AnnDataSet(adatas = sample2files,
                                        filename = outsdsfnm,
                                        add_key = 'sample')
sds.close()

# cost 50G to load
sds: sa2.AnnDataSet = sa2.read_dataset(filename = outsdsfnm,
                                       mode = 'r')
outadata_fnm = os.path.join(proj_dir, "18.snap2_peakcalling",
                            "out/scfilter/sa2pmat.union_pmat.adata.h5ad")
adata = sds.to_adata(copy_x = True, file = outadata_fnm)
sds.close()

new_obs_names = modify_obs_name(adata, obs_key = "sample")
adata.obs_names = new_obs_names
adata.close()

# * get AnnDataSet and AnnData for rdn peak
sa2pmatd2 = os.path.join(sa2pmatd, "union_pmat_rnd")
suffix = "_rnd_pmat.h5ad"
sample2files: List[Tuple[str, str]] = [
    (sample, os.path.join(sa2pmatd2, f"{sample}{suffix}")) for sample in samples
]
outsdsfnm = os.path.join(proj_dir, "18.snap2_peakcalling",
                         "out/scfilter/sa2pmat.pmat_rnd.h5ad")
if not os.path.exists(os.path.dirname(outsdsfnm)):
    os.makedirs(os.path.dirname(outsdsfnm), exist_ok = True)
    
sds: sa2.AnnDataSet = sa2.AnnDataSet(adatas = sample2files,
                                        filename = outsdsfnm,
                                        add_key = 'sample')
sds.close()

sds: sa2.AnnDataSet = sa2.read_dataset(filename = outsdsfnm,
                                       mode = 'r')
outadata_fnm = os.path.join(proj_dir, "18.snap2_peakcalling",
                            "out/scfilter/sa2pmat.pmat_rnd.adata.h5ad")
adata = sds.to_adata(copy_x = True, file = outadata_fnm)
sds.close()
new_obs_names = modify_obs_name(adata, obs_key = "sample")
adata.obs_names = new_obs_names
adata.close()

## make a copy of adata for rnd and union
## at /projects/ps-renlab/szu/

# * get fraction of cells for peaks
# * rnd peak
adata_rnd = sa2.read(os.path.join(proj_dir, "18.snap2_peakcalling",
                              "out/scfilter/sa2pmat.pmat_rnd.adata.h5ad"),
                 backed = "r")
adata_rnd = adata_rnd.to_memory()
groupby_rnd = atacmeta.loc[adata_rnd.obs_names, "pL4"]
pL4s = groupby_rnd.unique().tolist()

def get_peakfrac_rnd(group: str):
    t = adata_rnd.X[groupby_rnd == group, :]
    r = (t != 0).sum(axis = 0) / t.shape[0]
    return r
# fast, but each will cost 20G
with Pool(10) as p:
    peakfrac_rnd = p.map(get_peakfrac_rnd, pL4s)
peakfrac_rnd_mat = np.concatenate(peakfrac_rnd, axis = 0)
# to pandas
peakfrac_rnd_df = pd.DataFrame(peakfrac_rnd_mat, index = pL4s)
# save to pickle
peakfrac_rnd_df.to_pickle(
    os.path.join(proj_dir, "18.snap2_peakcalling", "out/scfilter",
                 "peakfrac_rnd.pkl"))

# * union peak
adata_union = sa2.read(os.path.join(proj_dir, "18.snap2_peakcalling",
                                    "out/scfilter/sa2pmat.union_pmat.adata.h5ad"),
                       backed = "r")
adata_union = adata_union.to_memory()

groupby_union = atacmeta.loc[adata_union.obs_names, "pL4"]
pL4s = groupby_union.unique().tolist()

def get_peakfrac_union(group: str):
    t = adata_union.X[groupby_union == group, :]
    r = (t != 0).sum(axis = 0) / t.shape[0]
    return r
# 2 * each cost 150G
# with additional one with 150G
with Pool(2) as p:
    peakfrac_union = p.map(get_peakfrac_union, pL4s)
peakfrac_union_mat = np.concatenate(peakfrac_union, axis = 0)
# to pandas
peakfrac_union_df = pd.DataFrame(peakfrac_union_mat, index = pL4s)
# set colnames
peaknms = adata_union.var_names.tolist()
peakfrac_union_df.columns = peaknms
# save to pickle
peakfrac_union_df.to_pickle(
    os.path.join(proj_dir, "18.snap2_peakcalling", "out/scfilter",
                "peakfrac_union.pkl"))

# * single-cell level pmat
# when pmat without any description, it means count
barcode2subclass: pd.Series = atacmeta[
    'subclass_label_v3'].apply(rename_allensubclass)

# load finalized peaks
final_peak_bed_fnm = os.path.join(
    proj_dir, "supple.07.peakcalling.allinone",
    "mba.whole.sa2.final.peak.srt.bed")
final_peaks: pd.DataFrame = read_bed4(final_peak_bed_fnm)


# * load AnnData with pmat of union peak
adata_union = sa2.read(
    os.path.join(proj_dir, "18.snap2_peakcalling",
                 "out/scfilter",
                 "sa2pmat.union_pmat.adata.h5ad"),
    backed = 'r'
)
# 120G memory fo snap if loaded to memory
# once adata loaded into memory,
# it becomes adata.AnnData, not SnapATAC2 AnnData
# so we can use functions from adata.AnnData
# snap = adata_union.to_memory()
# check if all the final peaks are there.
adata_var_mask = pd.Series(adata_union.var_names).isin(final_peaks.name)
adata_obs_group = barcode2subclass[adata_union.obs_names]
nds = 5000
seed = 0
def myds(data):
    if len(data) <= nds:
        return data
    return(data.sample(n = nds, replace = False, random_state = see
d))
adata_obs_group_ds = adata_obs_group.groupby(adata_obs_group).apply(myds)
adata_obs_ds = adata_obs_group_ds.index.get_level_values('barcode2')
outfnm_sa2ds = os.path.join(proj_dir, "18.snap2_peakcalling",
                            "out/scfilter",
                            "sa2pmat.final.peak.srt.ds5000.h5ad")
sa2_ds = adata_union.subset(obs_indices = adata_obs_ds,
                            var_indices = final_peaks.name,
                            out = outfnm_sa2ds)
sa2_ds.obs['subclassv3'] = adata_obs_group_ds

# to anndata's AnnData object
ann_ds: anndata.AnnData = sa2_ds.to_memory()
sa2_ds.close()
outfnm_ann2ds = os.path.join(proj_dir, "18.snap2_peakcalling",
                            "out/scfilter",
                            "ann.pmat.final.peak.srt.ds5000.h5ad")
ann_ds.write_h5ad(outfnm_ann2ds)

# * get subclass-level CPM
## about 35G in RAM
annds: anndata.AnnData = anndata.read_h5ad(
    os.path.join(proj_dir, "18.snap2_peakcalling",
                 "out/scfilter",
                 "ann.pmat.final.peak.srt.ds5000.h5ad"))
npeak = annds.shape[1]
ann_meta = annds.obs
grouped = ann_meta.groupby("subclassv3")
nsc = grouped.ngroups
cnt_pbysc = pd.DataFrame(
    np.zeros( (npeak, nsc) , dtype = np.float64),
    columns = list(grouped.groups.keys()),
    index = annds.var_names)

for group, idx in grouped.indices.items():
    cnt_pbysc[group] = np.ravel(
        annds.X[idx, ].sum(axis = 0, dtype = np.float64))

cnt_pbysc.to_csv(
    os.path.join(proj_dir, "18.snap2_peakcalling",
                 "out/scfilter",
                 "count_peakBysubclass.csv"),
    index = True,
    header = True
)



sum_pofsc = cnt_pbysc.sum(axis = 0)

values = cnt_pbysc.values
factors = sum_pofsc.values

cpm_pbysc = (values / factors ) * 10e6
cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
cpm_pbysc_df.index = cnt_pbysc.index
cpm_pbysc_df.to_csv(
    os.path.join(proj_dir, "18.snap2_peakcalling",
                 "out/scfilter",
                 "cpm_peakBysubclass.csv"),
    index = True,
    header = True
)





