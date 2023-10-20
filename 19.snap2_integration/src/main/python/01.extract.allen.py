import os
import sys
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from anndata import AnnData

# * configs
tscc_path: str = "/projects/ps-renlab2/szu/projects"
allen_dir: str = os.path.join(tscc_path, "shared_data/allen")
out_meta_dir: str = os.path.join(
    tscc_path, "CEMBA2/19.snap2_integration/src/main/resource")
# * work on annotation
# ** read annotation files
annot = pd.read_table(os.path.join(tscc_path, "shared_data/allen/AIT21_annotation.tsv"),
                      sep = "\t",
                      header = 0,
                      index_col = False)
annot.set_index('cluster_id', drop = False, inplace = True)

# ** map cluster/suptertype/subclass to nn or neuron
@dataclass
class AllenGroupLabel:
    cl: int
    id: int
    level: str
    label: str
    label_pre: int | None
    label_mid: str
    label_suf: int | None
    character: str
    def __hash__(self):
        return hash(self.id)


def match_group_label(label: str) -> Tuple[int|None, str, int|None]:
    match label:
        case [x, y, z]:
            return (int(x), y, int(z))
        case [x, y]:
            return (None, x, int(y))
        case [x]:
            return (None, x, None)
            
def init_AllenGroupLabel(cluster_id: int,
                label_level: str,
                annot: pd.DataFrame = annot, s1 = "_", s2 = " ") -> AllenGroupLabel:
    cl = annot.loc[cluster_id]["cl"]
    id = annot.loc[cluster_id][f"{label_level}_id"]
    label = annot.loc[cluster_id][f"{label_level}_label"]
    a, b, c = match_group_label(label.split(s1))
    
    clabel: AllenGroupLabel = AllenGroupLabel(
        cl = cl,
        id = id,
        level = label_level,
        label = label,
        label_pre = a,
        label_mid = b,
        label_suf = c,
        character = b.split(s2)[-1]
    )
    return(clabel)

def isGoodQuality(cl: AllenGroupLabel) -> bool:
    if cl.character == "LQ":
        return False
    return True

def init_List_of_AllenGroupLabel(label_level: str,
                                 annot: pd.DataFrame = annot,
                                 filter_LQ: bool = True,
                                 s1 = "_", s2 = " ") -> List[AllenGroupLabel]:
    r1 = map(lambda x: init_AllenGroupLabel(x, label_level, annot, s1,s2),
             annot.index)
    if filter_LQ:
        print("filter LQ.")
        r1 = filter(isGoodQuality, r1)
    return sorted(set(r1), key = lambda s: s.id)

def isnn(cl: AllenGroupLabel) -> bool:
    """non-neuronal cells and immmature neurons"""
    if cl.character in ["NN", "IMN"]:
        return True
    return False

cluster_labels = init_List_of_AllenGroupLabel("cluster")
cluster_labels_all = init_List_of_AllenGroupLabel("cluster", filter_LQ = False)
# supertype_labels = init_List_of_AllenGroupLabel("supertype")
# subclass_labels = init_List_of_AllenGroupLabel("subclass")
# class_labels = init_List_of_AllenGroupLabel("class")

# get nn and neuron related cl
nn_cls: set[int] = {i.cl for i in cluster_labels if isnn(i)}
neuron_cls: set[int] = {i.cl for i in cluster_labels if not isnn(i) }
lq_cls :set[int] = {i.cl for i in cluster_labels_all if not isGoodQuality(i)}

# get cl to cluster_id
cl2cid: Dict[int, int] = {
    row['cl'] : row['cluster_id'] for _, row in annot.iterrows()}

# * 10xv3 data
# load with memory-save mode
# 132G in RAM
# AnnData object
# with n_obs × n_vars = 2342082 × 32285 
#     obs: 'library_prep', 'gene.counts.0', 'doublet_score', 'roi',
#           'umi.counts', 'method', 'sex', 'external_donor_name',
#           'age', 'medical_conditions', 'cl'
#     var: 'gene_symbol', 'gene_identifier'
#     obsm: 'EMBED2', 'EMBED3'
#     layers: 'rawcount'
allen_10xv3 = ad.read_h5ad(
    filename = os.path.join(allen_dir, "AIT21_10Xv3.h5ad"),
    backed = 'r')

# * handle h5ad
def filter_obs(obs: pd.DataFrame, incls: set[str], sex = "M"):
    cond = (obs['sex'] == sex) & (obs['cl'].astype(int).isin(incls))
    return obs[cond]

def save_ann_with_filter(ann: AnnData,
                          filter_obs: pd.DataFrame,
                          outfnm: str,
                          colnm: str = "barcode") -> None:
    print(f"{filter_obs.shape[0]} barcodes in filter_obs.")
    r = ann[ann.obs[colnm].isin(filter_obs[colnm])].copy(outfnm)
    print(f"write filtered ann data to: {outfnm}")
    print(f"{r.shape[0]} barcodes in filtered ann.")
    return None

def check_obs(obs: pd.DataFrame,
              nn_cls: set[int] = nn_cls,
              neuron_cls: set[int] = neuron_cls,
              lq_cls: set[int] = lq_cls,
              sex: str = "M") -> Dict[str, int]:
    total_barcodes = obs.shape[0]
    male_barcodes = sum(obs['sex'] == sex)
    male_nn = sum( (obs['sex'] == sex) & (obs['cl'].astype(int).isin(nn_cls)) )
    male_neuron = sum( (obs['sex'] == sex) & (obs['cl'].astype(int).isin(neuron_cls)) )
    male_lq = sum( (obs['sex'] == sex) & (obs['cl'].astype(int).isin(lq_cls)) )
    r = {
        "all": total_barcodes,
        "male": male_barcodes,
        "male_nn": male_nn,
        "male_neuron": male_neuron,
        "male_lq": male_lq}
    a_sum = r["male_nn"] + r["male_neuron"] + r["male_lq"]
    if a_sum == r["male"]:
        print("number of barcodes in male matches.")
    else:
        print("number of barcodes in male does not match.")
    return r
    
    

obs_10xv3 = allen_10xv3.obs
obs_10xv3['barcode'] = obs_10xv3.index
obs_10xv3.to_csv(os.path.join(out_meta_dir, "allen_10xv3_obs.csv"), index = False)

# check to_csv, if cl is the real one or the faked categories index.
# obs_10xv3_csv = pd.read_csv(os.path.join(out_meta_dir, "allen_10xv3_obs.csv"),
#                             sep = ",",
#                             header = 0)
# cl_csv = obs_10xv3_csv["cl"]
# cl_obs = obs_10xv3["cl"].astype(int)
# a = [a1 == a2 for a1, a2 in zip(cl_csv, cl_obs)]
# so csv store the raw cl

nn_male_10xv3_obs = filter_obs(obs_10xv3, nn_cls)
neuron_male_10xv3_obs = filter_obs(obs_10xv3, neuron_cls)

check_10xv3 = check_obs(obs_10xv3, nn_cls, neuron_cls, lq_cls, sex = 'M')


save_ann_with_filter(allen_10xv3,
                     filter_obs = nn_male_10xv3_obs,
                     outfnm = f"{out_meta_dir}/nn_male_10xv3_ann.h5ad")
save_ann_with_filter(allen_10xv3,
                     filter_obs = neuron_male_10xv3_obs,
                     outfnm = f"{out_meta_dir}/neuron_male_10xv3_ann.h5ad")

# * 10xv2 data
allen_10xv2 = ad.read_h5ad(
    filename = os.path.join(allen_dir, "AIT21_10Xv2.h5ad"),
    backed = 'r'
)
obs_10xv2 = allen_10xv2.obs
obs_10xv2['barcode'] = obs_10xv2.index
obs_10xv2.to_csv(os.path.join(out_meta_dir, "allen_10xv2_obs.csv"), index = False)
nn_male_10xv2_obs = filter_obs(obs_10xv2, nn_cls)
neuron_male_10xv2_obs = filter_obs(obs_10xv2, neuron_cls)

check_10xv2 = check_obs(obs_10xv2, nn_cls, neuron_cls, lq_cls, sex = 'M')

save_ann_with_filter(allen_10xv2,
                     filter_obs = nn_male_10xv2_obs,
                     outfnm = f"{out_meta_dir}/nn_male_10xv2_ann.h5ad")
save_ann_with_filter(allen_10xv2,
                     filter_obs = neuron_male_10xv2_obs,
                     outfnm = f"{out_meta_dir}/neuron_male_10xv2_ann.h5ad")


# * snRNA data
allen_multiome = ad.read_h5ad(
    filename = os.path.join(allen_dir, "AIT21_multiome.h5ad"),
    backed = 'r'
)
obs_multiome = allen_multiome.obs
obs_multiome['barcode'] = obs_multiome.index
obs_multiome.to_csv(os.path.join(out_meta_dir, "allen_multiome_obs.csv"), index = False)
nn_male_multiome_obs = filter_obs(obs_multiome, nn_cls)
neuron_male_multiome_obs = filter_obs(obs_multiome, neuron_cls)

check_multiome = check_obs(obs_multiome, nn_cls, neuron_cls, lq_cls, sex = 'M')

save_ann_with_filter(allen_multiome,
                     filter_obs = nn_male_multiome_obs,
                     outfnm = f"{out_meta_dir}/nn_male_multiome_ann.h5ad")
save_ann_with_filter(allen_multiome,
                     filter_obs = neuron_male_multiome_obs,
                     outfnm = f"{out_meta_dir}/neuron_male_multiome_ann.h5ad")




