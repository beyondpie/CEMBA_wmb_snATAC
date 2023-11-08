import pickle
import os
import sys
import numpy as np
import pandas as pd

with open(f"../result/clustering_sum_L1/sa2_clustering_L0_0.pkl", 'rb') as f:
    L1_sum = pickle.load(f)

with open(f"../result/clustering_sum_L1/barcodes.txt", 'r') as f:
    barcodes = [l.strip() for l in f.readlines()]

r = 0.4
which_col = L1_sum['leiden_r'] == r
leiden = L1_sum['leiden'][:, which_col]

with open(f"../result/clustering_sum_L1/sa2_L1_r0.4_barcodes2id.csv", 'w') as f:
    f.writelines("barcode,L1\n")
    f.writelines([f"{b},{i}\n" for b, i in zip(barcodes,leiden[:,0])])

# copy this file to 17.snapatac2/resource/sa2_dlt2_barcode2id.csv

# save umap
umap = L1_sum['umap']
barcode2umap = pd.DataFrame(data = {"barcode" : barcodes,
                                    "UMAP1" : umap[:,0],
                                    "UMAP2" : umap[:,1]})
barcode2umap.to_csv("../result/clustering_sum_L1/L1_UMAP.csv",
                    header = True, index = False)
