"""
This script is used to record how to read some supplementary data,
which are published in Nature already.
"""
from typing import List
import os
import re
import numpy as np
import pandas as pd


# * Download Supplementary Data
# You can find the link under the paper website:
# https://www.nature.com/articles/s41586-023-06824-9
# Here I put it under my Download directory

# * get subclass-specific cCREs
# load L4-level cCRE first.
L4toCRE_fnm = os.path.join(
    "/Users/szu/Downloads",
    "2023-03-05453B-s3/SI Tables",
    "SI Table 6 Cell-subtype specific cCREs.csv",
)
with open(L4toCRE_fnm, "r") as f:
    CREs: List[str] = f.readline().strip().split(",")

hasCRE: np.ndarray = np.loadtxt(
    fname=L4toCRE_fnm,
    delimiter=",",
    dtype=np.uint8,
    skiprows=1,
    usecols=range(1, len(CREs) + 1),
)
L4names: np.ndarray = np.loadtxt(
    fname=L4toCRE_fnm, delimiter=",", dtype=str, skiprows=1, usecols=0
)
L4toCRE: pd.DataFrame = pd.DataFrame(data=hasCRE, index=L4names, columns=CREs)
pL4: List[str] = L4names.tolist()

# L4 to subclass
L4meta_fnm = os.path.join(
    "/Users/szu/Downloads",
    "2023-03-05453B-s3/SI Tables",
    "SI Table 4 L4-level annotation based on the integration analysis with scRNA-seq data including major region information.xlsx",
)

L4meta: pd.DataFrame = pd.read_excel(io=L4meta_fnm, header=0, sheet_name=0)
L4meta.set_index(keys="ATAC L4 cluster", drop=False, inplace=True)


# pL4 to subclass
def map_pL4_to_sc(apL4: str) -> str:
    t = apL4.replace("p", "")
    tt: List[str] = t.split(sep="_")
    if len(tt) < 2:
        return L4meta.loc[t]["Subclass label"]
    l3: str = re.sub("-\d+$", "", tt[0])
    l4s: List[str] = [f"{l3}-{i}" for i in tt[1:]]
    l4s.append(tt[0])
    scs: List[str] = L4meta.loc[l4s]["Subclass label"].unique().tolist()
    if len(scs) > 1:
        print(f"{apL4} has multiple subclasses.")
    return scs[0]


scs: List[str] = [map_pL4_to_sc(k) for k in pL4]

# subclass to CREs
L4toCRE.insert(0, column="sc", value=scs)
sctoCRE = L4toCRE.groupby(by="sc").apply(
    lambda x: x[CREs].sum(axis = 0))

# output subclass' CREs
outdir = os.path.join("/Users/szu/git-recipes",
                      "mouseBrainAtlas/cCRE_heatmap/data",
                      "subclass2CRE")
def outscCRE(sc: str) -> None:
    scnm = sc.replace('/', '-').replace(' ', '_')
    outfnm = os.path.join(outdir, f"{scnm}.cCREs.txt")
    pd.Series(sctoCRE.columns[sctoCRE.loc[sc] > 0]).to_csv(
        outfnm, sep = "\t", header = False, index = False
    )

for sc in sctoCRE.index.tolist():
    print(f"{sc}")
    outscCRE(sc)
