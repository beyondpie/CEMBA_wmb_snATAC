import os
import sys
from typing import List
from pathlib import Path
import itertools

import numpy as np
import matplotlib.pyplot as plt
import snapatac2 as sa2

with open("../meta/mba.whole.sample.lst", 'r') as f:
    samples = [ l.strip() for l in f.readlines()]

# * summarize doublet removal under bmat
bmat_dlt_dir = "../result/barcode2dltprob"
dlt_prob_threshold = 0.5
sample2dlt = {}
for s in samples:
    with open(f"{bmat_dlt_dir}/{s}.txt", 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        barcode2dltprob = [
            (l.split(',')[0], float(l.split(',')[1])) for l in lines ]
    sample2dlt[s] = barcode2dltprob
barcodes_all_bmat = list(itertools.chain.from_iterable(
    [sample2dlt[s] for s in samples] ))  

sample2barcodes = {}
for s in samples:
    sample2barcodes[s] = [
        v[0] for v in list(filter(lambda x: x[1] <= dlt_prob_threshold,
                                  sample2dlt[s]))]
# 2355842
nbarcodes_after_dlt = sum(
    [len(sample2barcodes[s]) for s in sample2barcodes.keys()])

# * load results from doublet removal under gmat
gmat_dlt_dir = "../../00.data.preprocess/snapatac2_pp_out/pp_stat"
sample2barcodes_gmat = {}
for s in samples:
    with open(f"{gmat_dlt_dir}/{s}.qc.dlt.barcodes", 'r') as f:
        sample2barcodes_gmat[s] = [l.strip() for l in f.readlines()]

# 2361710
nbarcodes_gmat = sum(
    [len(sample2barcodes_gmat[s]) for s in sample2barcodes_gmat.keys()])

# ** joint between barcodes under bmat and gmat
barcodes_bmat = list(itertools.chain.from_iterable(
    [sample2barcodes[s] for s in sample2barcodes.keys()]))
barcodes_gmat = list(itertools.chain.from_iterable(
    [sample2barcodes_gmat[s] for s in sample2barcodes_gmat.keys()]
))

# 2326359
barcodes_both = list(
    set(barcodes_bmat).intersection(set(barcodes_gmat)))

len(barcodes_both) / len(barcodes_bmat) # 98.75%

# * compare with SnapATAC
with open("../../supple.02.QC/sa1.qc.dlt.barcodes", 'r') as f:
    barcodes_sa1 = [l.strip() for l in f.readlines()]
# 2204291
barcodes_sa1_bmat = list(
    set(barcodes_sa1).intersection(set(barcodes_bmat))
)

# * draw dlt rate
dl2r = {}
for s in samples:
    dl2r[s] = 1 - len(sample2barcodes[s]) / len(sample2dlt[s])

with open("../../supple.02.QC/sample2biorep.csv", 'r') as f:
    lines = [l.strip() for l in f.readlines()]
    s2rep = {l.split(',')[0]: l.split(',')[1] for l in lines}

dlt_early = [dl2r[s] for s in samples if s2rep[s] == 'early']
dlt_later = [dl2r[s] for s in samples if s2rep[s] == 'later']
fig = plt.figure()
plt.boxplot([dlt_early, dlt_later])
plt.show()

s_early = [s for s in samples if s2rep[s] == 'early']
s_later = [s for s in samples if s2rep[s] == 'later']
