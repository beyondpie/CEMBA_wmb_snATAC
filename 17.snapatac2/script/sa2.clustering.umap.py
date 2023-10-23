import os
import sys
import math
import pickle
from typing import Dict, List
from dataclasses import dataclass, field
import numpy as np
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
from leiden import LeidenSum, ScatterPlot
from leiden import draw_umap
from leiden import init_LeidenSum_from_file
from colors import SnapATACPalette

    
if __name__ == '__main__':
    from_dir: str = sys.argv[1]
    cll: str = sys.argv[2]
    outdir: str = sys.argv[3]
    nsample: int = int(sys.argv[4])
    # this will be used in parallel
    cid: str = sys.argv[5]
    ls: LeidenSum = init_LeidenSum_from_file(
        from_dir = from_dir,
        cll = cll,
        cid = cid)
    print(f"draw UMAP for {cll}: {cid}")
    scf = ScatterPlot(nsample = nsample)
    ## FIXME: put colors into scf
    draw_umap(t = ls, scf = scf, outdir = outdir, colors = list(set(SnapATACPalette)))

