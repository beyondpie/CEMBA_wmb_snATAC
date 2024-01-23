from typing import List
import snapatac2 as sa2


def modify_obs_name(sds: sa2.AnnData | sa2.AnnDataSet,
                    obs_key = "sample") -> List[str]:
    obs_names: List[str] = [f"{i}.{j}"
              for i, j in zip(
                      sds.obs[obs_key].to_list(), sds.obs_names)]
    return obs_names

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
