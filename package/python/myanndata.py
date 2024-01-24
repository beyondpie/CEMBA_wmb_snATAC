import warnings
warnings.filterwarnings('default', category=DeprecationWarning)
def script_is_deprecated():
    warnings.warn("myanndata.py is deprecated. Functions move to cembav2env.sa2.py.")
script_is_deprecated()


from typing import List, Dict
import numpy as np
import pandas as pd
import anndata as ad

def grouped_obs_mean(adata: ad.AnnData,
                     group_meta: pd.DataFrame,
                     group_key: str,
                     layer: str|None=None,
                     gene_symbols: List[str]|None =None) -> pd.DataFrame:
    """
    Ref: https://github.com/scverse/scanpy/issues/181
    Add checks before running the content.
    """
    # check
    if (layer is not None) and (layer not in adata.layers.keys()):
        raise KeyError(f"{layer} not found in adata layer.")
    # group_meta should have two columns: barcode, group_key
    if ("barcode" not in group_meta.columns) or (group_key not in group_meta.columns):
        raise KeyError(f"Either barcode or {group_key} not in group_meta.")
    # and group_meta index by barcode
    if (adata.shape[0] != group_meta.shape[0]):
        raise RuntimeError("adata and group_meta have diff # of cells.")
    # TODO: check the order of group_meta is the as adata's barcode
    
    # make sure the anndata and group_meta have the same order of barcodes
    group_meta.set_index('barcode', drop = False)
    group_meta = group_meta.loc[adata.obs_names]

    idy = (
        [True]
        if gene_symbols is None
        else adata.var_names.isin(gene_symbols)
    )
    if sum(idy) == 0:
        raise KeyError("No gene_symbols found in adata.")
    match gene_symbols:
        case None:
            subset_adata = adata
        case _:
            print(f"{sum(idy)} of {len(gene_symbols)}",
                  "genes found in adata.")
            subset_adata = adata[:, idy]
    match layer:
        case None:
            X = subset_adata.X
        case _:
            print(f"use {layer} for X.")
            X = subset_adata.layers[layer]
    grouped = group_meta.groupby(group_key)
    print(f"find {grouped.ngroups} under group {group_key}")
    out = pd.DataFrame(
        np.zeros((subset_adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=subset_adata.var_names
    )
    for group, idx in grouped.indices.items():
        out[group] = np.ravel(X[
idx, :].mean(axis=0, dtype=np.float64))
    return out

def get_barcode2group(obs: pd.DataFrame,
                      annot: pd.DataFrame,
                      fromcol: str = "cl",
                      tocol: List[str] = ["subclass_id"],
                      index_col: str = "barcode"):
    """Generate by GPT-4"""
    # Create a DataFrame with the mapping
    f2t = annot.loc[:, [fromcol] + tocol].copy()

    from_obs = obs[fromcol]
    if from_obs.dtype.name == 'category':
        print(f"{fromcol} is category, change to int")
        obs[fromcol] = from_obs.astype(int)

    # Merge the mapping DataFrame with the original obs DataFrame
    result = obs.merge(f2t, how='left', on=fromcol)
    result.set_index(index_col, drop = False, inplace = True)
    return result
