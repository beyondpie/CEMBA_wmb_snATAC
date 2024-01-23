from inspect import getsource
import logging
from typing import List, Dict
import warnings
import numpy as np
import pandas as pd

def printfn(fn) -> None:
    """print fn source code"""
    content = getsource(fn)
    print(content)

def concat_allen_name(allen_l2):
    """allen_l2 is usually a pandas Series.
    """
    import re
    r = re.sub("_$", "", allen_l2.replace("/", "-").replace(" ", "_"))
    return(r)

def set_file_logger(fnm:str,
                    fmode:str = 'a',
                    name:str = 'sa2_pp',
                    log_level: int = logging.DEBUG) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(filename = fnm,
                             mode = fmode)
    fm = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fm)
    logger.addHandler(fh)
    return logger

def head(obj, n = 20):
    if hasattr(obj, "head"):
        return obj.head(n)
    if type(obj) in [list, tuple]: 
        return(obj[:n])
    if hasattr(obj, "keys"):
    # if type(obj) is dict:
        keys = list(obj.keys())[:n]
        return({k:obj[k] for k in keys})
    if hasattr(obj, "__iter__"):
        tmp = 0
        a = iter(obj)
        r = []
        try:
            while tmp < n:
                tmp += 1
                r.append(next(a))
        finally:
            return r
    else:
        print(f"no support for {type(obj)}.")
        return None

def clean_AnnDataSet(d: str):
    import os
    os.remove(f"{d}/_dataset.h5ads")
    import shutil
    shutil.rmtree(f"{d}/anndatas", ignore_errors = True)
    os.removedirs(d)

# deprecated
# use to_adata directly for subsetting data.
def get_subset_from_AnnDataSet(adata, outf,
                               obs_index = None, var_index = None,
                               logger = None):
    import os
    tmp_dir, _ = os.path.splitext(outf)
    os.makedirs(tmp_dir, exist_ok = True)
    adata_subset = adata.subset(
        obs_indices = obs_index,
        var_indices = var_index,
        out = tmp_dir
    )
    adata_subset = adata_subset[0]
    if logger is not None:
        logger.info(f"create subset from {tmp_dir} to {outf}")
    r = adata_subset.to_adata(
        copy_x = True,
        file = outf
    )
    adata_subset.close()
    # remove adata_subset
    os.remove(f"{tmp_dir}/_dataset.h5ads")
    import shutil
    shutil.rmtree(f"{tmp_dir}/anndatas", ignore_errors = True)
    os.removedirs(tmp_dir)
    if logger is not None:
        logger.info(f"clean tmp {tmp_dir}.")
    return r


def modify_obs_name(sds, obs_key = "sample") -> List[str]:
    obs_names: List[str] = [f"{i}.{j}"
              for i, j in zip(
                      sds.obs[obs_key].to_list(), sds.obs_names)]
    return obs_names
    
def deprecated(func):
    """
    Usage:
    @deprecated
    def old_function():
        print("This is the old function.")
    """
    def wrapper(*args, **kwargs):
        warnings.warn(f"{func.__name__} is deprecated.",
                      DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    return wrapper


def complete_values(key2value: Dict[str, str],
                       all_keys: List[str],
                       default_value: str = "NA") -> List[str]:
    keys_in_pool = set(key2value.keys())
    is_in_pool = [b in keys_in_pool for b in all_keys]
    if all(is_in_pool):
        print("all_keys in key2value.")
    else:
        keys_left = [k for i, k in enumerate(all_keys) if not is_in_pool[i]]
        print(f"{len(keys_left)} not in key2value.")
        for k in keys_left:
            key2value[k] = default_value
    return [key2value[k] for k in all_keys]


def shuffle_column(df: pd.DataFrame,
                  bycol:str = 'rep',
                  random_seed:int = 0) -> pd.Series:
    np.random.seed(random_seed)
    shuffled_index = np.random.permutation(df.index, )
    return df.loc[shuffled_index, bycol]
