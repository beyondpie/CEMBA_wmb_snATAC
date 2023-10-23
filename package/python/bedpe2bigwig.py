import gzip
import os
import pandas as pd
from typing import Dict, List

def generate_bedpe_worker(meta: str,
                          sample: str,
                          infnm: str,
                          outdir:str,) -> bool:
    """Partition barcodes from a sample bedpe to different goups.

    Outfile format:
    fnm = f"{outdir}/{sample}.{g}.bedpe.gz" 
    """
    # * read meta data
    if not infnm.endswith('.bedpe.gz'):
        raise FileNotFoundError(
            f"{infnm} should be end with .bedpe.gz")
    r: pd.DataFrame = pd.read_csv(meta, sep = "\t", header = 0)
    barcode = pd.Series(
        data = r.iloc[:, 2].values,
        index = r.iloc[:, 0])
    cell2group :Dict[str, str] = barcode.to_dict()
    # * set up outfiles
    s = r.iloc[:, [1,2]].groupby(by = r.iloc[:,1]).apply(
        lambda x: list(pd.unique(x.iloc[:,1])))
    o2n :Dict[str, List[str]] = s.to_dict()
    ng2f :Dict[str, gzip.GzipFile] = dict()

    if not os.path.exists(outdir):
        print(f"{outdir} does not exist.")
        os.makedirs(name = outdir, exist_ok = True)
    for g in o2n[sample]:
        fnm = f"{outdir}/{sample}.{g}.bedpe.gz"
        if os.path.exists(fnm):
            print(f"{fnm} exists and remove it.")
            os.remove(fnm)
        ng2f[g] = gzip.open(f"{fnm}", mode = "wb")

    # * main
    with gzip.open(infnm, "rb") as f:
        for l in f:
            barcode: str = l.decode().split("\t")[6]
            if barcode in cell2group:
                ng:str = cell2group[barcode]
                ## debug
                # print(f"{barcode} -> {ng}")
                ng2f[ng].write(l)
    # close all the output files
    for g in o2n[sample]:
        if not ng2f[g].closed:
            ng2f[g].close()
    return True
