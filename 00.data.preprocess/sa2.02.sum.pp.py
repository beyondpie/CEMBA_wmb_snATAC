import snapatac2 as sa2
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


sample_file = "mba.whole.sample.lst"
# * get all the barcodes after qc + doublets
with open(sample_file, 'r') as f:
    samples = [l.strip() for l in f.readlines()]

# input_dir = "snapatac2_pp_out/bmat"
# out_dir = "snapatac2_pp_out/pp_stat"
# if not os.path.exists(out_dir):
#     os.makedirs(out_dir, exist_ok = True)

# for sample in samples:
#     print(f"Sample: {sample}")
#     sd = sa2.read(f"{input_dir}/{sample}.sa2.bmat.hdf5")
#     barcodes = [f"{sample}.{b}" for b in sd.obs_names]
#     out_file = f"{out_dir}/{sample}.qc.dlt.barcodes"
#     with open(out_file, 'w') as f:
#         f.write('\n'.join(barcodes))
#         print(f"{sample} barcodes to file.")

# * barcodes after QC
# input_dir = "snapatac2_pp_out/qc"
# out_dir = "snapatac2_pp_out/qc_stat"
# if not os.path.exists(out_dir):
#     os.makedirs(out_dir, exist_ok = True)
# for sample in samples:
#     print(f"Sample: {sample}")
#     sd = sa2.read(f"{input_dir}/{sample}.sa2.qc.hdf5")
#     barcodes = [f"{sample}.{b}" for b in sd.obs_names]
#     out_file = f"{out_dir}/{sample}.qc.barcodes"
#     with open(out_file, 'w') as f:
#         f.write('\n'.join(barcodes))
#         print(f"{sample} barcodes to file.")

# * barcodes raw
input_dir = "snapatac2_pp_out/raw_fragment"
out_dir = "snapatac2_pp_out/raw_barcode"
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok = True)
for sample in samples:
    print(f"Sample: {sample}")
    sd = sa2.read(f"{input_dir}/{sample}.sa2.raw.hdf5")
    barcodes = [f"{sample}.{b}" for b in sd.obs_names]
    out_file = f"{out_dir}/{sample}.raw.barcodes"
    with open(out_file, 'w') as f:
        f.write('\n'.join(barcodes))
        print(f"{sample} barcodes to file.")

s2b_raw = {}
for s in samples:
    with open(f"snapatac2_pp_out/raw_barcode/{s}.raw.barcodes") as f:
        barcodes = [l.strip() for l in f.readlines()]
        s2b_raw[s] = set(barcodes)
n_barcode_raw = [len(s2b_raw[s]) for s in samples]


# * total number of barcodes afte pp.
# sample2barcodes = {}
# for s in samples:
#     with open(f"{out_dir}/{s}.qc.dlt.barcodes", 'r') as f:
#         barcodes = [l.strip() for l in f.readlines()]
#         sample2barcodes[s] = barcodes

# n_barcodes = [len(sample2barcodes[s])
#               for s in sample2barcodes.keys()]

# s2b_qc = {}
# for s in samples:
#     with open(f"snapatac2_pp_out/qc_stat/{s}.qc.barcodes") as f:
#         barcodes = [l.strip() for l in f.readlines()]
#         s2b_qc[s] = barcodes
n_barcode_qc = [len(s2b_qc[s]) for s in samples]

dl2r = {}
for s in samples:
    dl2r[s] = 1 - len(sample2barcodes[s]) / len(s2b_qc[s])

# barcodes_sa2 = []
# for s in samples:
#     barcodes_sa2.extend(sample2barcodes[s])
# # In total, we have 2361710 barcodes left.

# # * compare with SnapATAC
# with open("../supple.02.QC/sa1.qc.dlt.barcodes", 'r') as f:
#     barcodes_sa1 = [l.strip() for l in f.readlines()]
# barcodes_sa1_sa2 = list(set(barcodes_sa1) & set(barcodes_sa2))

# * dlt rate
dl2r = {}
for s in samples:
    dl2r[s] = 1 - len(sample2barcodes[s]) / len(s2b_qc[s])

dlr = [dl2r[s] for s in samples]

fig = plt.figure()
plt.boxplot(dlr)
plt.show()

# * read sampel2 biorep
with open("../supple.02.QC/sample2biorep.csv", 'r') as f:
    lines = [l.strip() for l in f.readlines()]
    s2rep = {l.split(',')[0]: l.split(',')[1] for l in lines}

dlt_early = [dl2r[s] for s in samples if s2rep[s] == 'early']
dlt_later = [dl2r[s] for s in samples if s2rep[s] == 'later']
fig = plt.figure()
plt.boxplot([dlt_early, dlt_later])
plt.show()

s_early = [s for s in samples if s2rep[s] == 'early']
s_later = [s for s in samples if s2rep[s] == 'later']


