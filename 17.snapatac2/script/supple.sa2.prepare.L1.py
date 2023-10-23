import os
import sys
from typing import List
from pathlib import Path
import itertools

import numpy as np
import snapatac2 as sa2

# get barcodes after qc and doublet removal
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

barcodes = []
for s in samples:
    barcodes.extend(sample2barcodes[s])

with open("../resource/barcode2id_L0.csv", 'w') as f:
    f.writelines("barcode,L0\n")
    f.writelines('\n'.join([f"{v},0" for v in barcodes]))

