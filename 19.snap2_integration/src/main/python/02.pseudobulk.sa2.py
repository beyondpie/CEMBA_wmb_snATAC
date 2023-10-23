import os
import sys

import pandas as pd
import anndata as ad
from pyprojroot import here
import re

proj_root = str(here())
sys.path.insert(0, os.path.join(proj_root, "package/python"))
from myanndata import grouped_obs_mean, get_barcode2group
from utils import set_file_logger

# * configs
resource_dir = f"{proj_root}/19.snap2_integration/src/main/resource"
logger = set_file_logger(fnm=f"{resource_dir}/avgexp_sa2.log", name="avgexp")

# * normalize raw gmat-based snap data 
## load snap data to ad.AnnData
# snap = sa2.read(
#     f"{proj_root}/17.snapatac2/sa2_sa2default_gmat/sa2default_gmat_merged.h5ad",
#     backed = 'r'
# )
# # 95G
# # takes about 20 minutes or more
# mysnap: ad.AnnData = snap.to_memory()
# snap.close()
# mysnap.obs['L2'] = b2g.loc[mysnap.obs_names]['L2']
# ## perform log CPM normlaization
# # this cost half an hour or even longer
# mysnap.layers['sa2default_raw'] = mysnap.X.copy()
# # this cost half an hour, and 250g mem
# sc.pp.normalize_total(mysnap, target_sum = 1e4)
# # this is fast
# sc.pp.log1p(mysnap)
# # ~195G in disk
# # takes about 10 minutes
# mysnap.write(f"{resource_dir}/sa2_gmat/sa2_sa2default_ann_with_CPMlognorm.h5ad")
# mysnap_tf = mysnap[:, TFs]
# # about 4G
# mysnap_tf.write(
#     f"{resource_dir}/sa2_gmat/sa2_sa2default_tf_ann_with_CPMlognorm.h5ad")
# # takes the similiar time and mem like above and
# # with error since need to write to tmp file
# # snap_tf = snap.subset(
# #     var_indices = pd.Series(snap.var_names).isin(TFs), out = None)
# # snap.close()

# ## load barcoede2group file.
# b2g = pd.read_csv(
#     f"{proj_root}/17.snapatac2/resource/sa2_dlt2_L2_barcode2id.csv",
#     header = 0, sep = ",")
# b2g.set_index('barcode', inplace = True, drop = True)

# mysnap_tf = ad.read_h5ad(
#     f"{resource_dir}/sa2_gmat/sa2_sa2default_tf_ann_with_CPMlognorm.h5ad")
# avg_tf_snap = grouped_obs_mean(mysnap_tf, group_key = "L2")
# avg_tf_snap.to_csv(
#     os.path.join(f"{resource_dir}/sa2_gmat",
#                  "snap_all_tf_by_scid_avg_lognorm.csv"))

# * load snap after normalization.
snap_file = os.path.join(resource_dir,
                              "sa2_gmat",
                              "sa2_sa2default_ann_with_CPMlognorm.h5ad")
logger.info(f"Load {snap_file} after lognorm into memory.")
mysnap: ad.AnnData = ad.read_h5ad(snap_file, backed = None)
logger.info("Done: to memory as anndata.AnnData instance.")

barcodes = mysnap.obs.index
# * load barcode2group file.
b2g_file = os.path.join(proj_root, "17.snapatac2", "post_script",
                        "sa2.barcode2L3.csv")
logger.info(f"Load b2g_file: {b2g_file}")
b2g:pd.DataFrame = pd.read_csv(b2g_file, sep = ",", header = 0)
b2g["L2"] = b2g["L3"].map(lambda x : re.sub(r"-[0-9]+$", "", x))
b2g["L1"] = b2g["L3"].map(lambda x: re.sub(r"-[0-9]+-[0-9]+$", "", x))
b2g.set_index("barcode", drop = False, inplace = True)
# IMPORTANT: make sure b2g have the same order as mysnap
b2g = b2g.loc[barcodes]

# check
# b2g["L1"].value_counts()
# b2g["L2"].value_counts()
# b2g["L3"].value_counts()

# get avg exp
groups = ["L1", "L2", "L3"]
out_dir = f"{resource_dir}/sa2_gmat"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
for g in groups:
    logger.info(f"Calculate avg for group: {g}.")
    avgexp = grouped_obs_mean(adata = mysnap,
                              group_meta= b2g,
                              group_key = g)
    avgexp.to_csv(
        os.path.join(out_dir, f"sa2_sa2default_{g}.avgexp.csv"))
    
