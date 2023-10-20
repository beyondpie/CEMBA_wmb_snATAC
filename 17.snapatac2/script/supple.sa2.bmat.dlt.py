import os
import sys
from typing import List
from pathlib import Path

import numpy as np
import snapatac2 as sa2

# test
# work_dir = "/Users/szu/git-recipes/mouseBrainAtlas/CEMBA2"
# sample = "CEMBA190718_8F"
# qc_dlt_dir = f"{work_dir}/17.snapatac2/result/qc_bmat-dlt/"
# out_dir = f"{work_dir}/17.snapatac2/result/qc_bmat-dlt_barcodes"

work_dir = "/oasis/tscc/scratch/szu/projects/CEMBA2/"
qc_dlt_dir = f"{work_dir}/17.snapatac2/sa2_qc_dlt/qc_dlt"
sample = sys.argv[1]
out_dir = f"{work_dir}/17.snapatac2/sa2_qc_dlt/barcode2dltprob"

os.makedirs(out_dir, exist_ok = True)
dlt_prob_threshold = 0.5

outf = f"{out_dir}/{sample}.txt"
snap = sa2.read(Path(f"{qc_dlt_dir}/{sample}_qc_dlt.h5ad"), backed = 'r')
barcodes: np.ndarray = np.array([f"{sample}.{k}" for k in snap.obs_names])
dlt_probs = snap.obs['doublet_probability'].to_numpy().tolist()
# barcodes_filtered:List[str] = barcodes[dlt_probs <= dlt_prob_threshold].tolist()
with open(outf, 'w') as f:
    f.writelines('\n'.join([f"{b},{s}" for b, s in zip(barcodes, dlt_probs)]))
snap.close()











