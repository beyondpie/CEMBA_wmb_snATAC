import os
import sys
import numpy as np
import pandas as pd
import snapatac2 as sa2

clevel = "L1"
out_dir = f"sa2_dlt2_{clevel}_subsets"
os.makedirs(out_dir, exist_ok = False)
# load snap data into memory:
# - it consumes about 200G and loading process takes about 10 minutes.
# - tscc has slow IO, so use this to make subset faster.
# - this loading format will limit to read this file for another loading.
#   - not sure if: read by 'r', then use to_memory will better, since this
#     will create a new object based on the description
# FIXME: after this loading, no subset attribute.


sds_all = sa2.read("resource/merge_cemba_all.h5ad", None)

barcode2id: pd.DataFrame = pd.read_csv("resource/sa2_dlt2_barcode2id.csv")

cid = 0
snap_file = f"{out_dir}/sa2_dlt2_{clevel}_{cid}.h5ad"
sub_barcodes = barcode2id[barcode2id[clevel] == int(cid)]['barcode']
all_barcodes = sds_all.obs_names
a = set(sub_barcodes)
is_in_sub = np.array([b in a for b in all_barcodes])
print(f"Find {is_in_sub.sum()} barcodes in cemba.")
sds = sds_all.subset(
    obs_indices = is_in_sub,
    out = snap_file)





