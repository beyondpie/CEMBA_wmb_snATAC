import os
import sys
import traceback
from typing import Dict, Tuple
from pathlib import Path
import pickle
import shutil

import numpy as np
import pandas as pd

from numba.core.errors import NumbaDeprecationWarning 
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import snapatac2 as sa2
import pyprojroot

code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils  # noqa: E402
from leiden import cemba  # pyright: ignore # noqa: E401, E402, F401
from leiden import cal_silhouette  # pyright: ignore # noqa: E402, F401
from leiden import umap  # pyright: ignore # noqa: E402

# * log
logger = utils.set_file_logger(  # pyright: ignore
    fnm=snakemake.log[0], name="leiden.united"  # pyright: ignore # noqa: F821
)
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

# * meta
snap_file = snakemake.output["snap_file"][0] # pyright: ignore # noqa
# embed
blacklist_file = snakemake.input["blacklist"]  # pyright: ignore # noqa
barcode2id_file = snakemake.input["barcode2id"]  # pyright: ignore # noqa
# cemba_anndataset_file = snakemake.input["cemba_anndataset"]  # pyright: ignore # noqa
cemba_anndata_file = snakemake.input["cemba_anndata"]  # pyright: ignore # noqa

embed_params: Dict = snakemake.params["embed"]  # pyright: ignore # noqa
embed_nm: str = embed_params["name"]
# knn
knn_params: Dict = snakemake.params["knn"]  # pyright: ignore # noqa
knn_nm: str = knn_params["name"]
# leiden
leiden_params: Dict = snakemake.params["leiden"]  # pyright: ignore # noqa
leiden_nm: str = leiden_params["name"]
umap_params: Dict = snakemake.params["umap"]  # pyright: ignore # noqa

# clustering info
cid = snakemake.params["cluster_id"]  # pyright: ignore # noqa
clevel = snakemake.params["clevel"]  # pyright: ignore # noqa

# out summary file
out_sum_file = f"{os.path.dirname(snap_file)}/sa2_clustering_{clevel}_{cid}.pkl"
clustering_sum: Dict = {
    "embed_params": embed_params,
    "knn_params": knn_params,
    "leiden_params": leiden_params,
    "umap_params": umap_params
}

# * main
# ** load cemba anndataset
logger.info(f"Loading cemba all AnnData: {cemba_anndata_file}.")
# use read only mode for reading the same data simutaneously
# sds_all = sa2.read_dataset(Path(cemba_anndataset_file), mode = 'r')
sds_all = sa2.read(Path(cemba_anndata_file), 'r')

# ** subset sds_all
logger.info(f"Subset above AnnData: {clevel}_{cid}.")
barcode2id: pd.DataFrame = pd.read_csv(barcode2id_file,
                                       header=0, index_col=False,
                                       dtype = str)
# cid here is str
sub_barcodes = barcode2id[barcode2id[clevel] == cid]['barcode']
logger.info(f"Under {clevel}_{cid}, find {sub_barcodes.size} barcodes.")
# cost lots of time when sub_barcodes and all_barcodes are large.
# is_in_sub = np.isin(all_barcodes, sub_barcodes, assume_unique = True)
# NOTE: alternate np.isin
# https://github.com/numpy/numpy/issues/14997
all_barcodes = sds_all.obs_names
a = set(sub_barcodes)
is_in_sub = np.array([b in a for b in all_barcodes])
logger.info(f"Find {is_in_sub.sum()} barcodes in cemba.")

# sds = sds_all.to_adata(
#     obs_indices = is_in_sub,
#     copy_x = True,
#     file = snap_file
# )
sds = sds_all.subset(
    obs_indices = is_in_sub,
    out = snap_file
)
logger.info(f"Use to_data to subet the sds_all to {snap_file}.")

# ** clustering procedures
logger.info("Select features.")
sa2.pp.select_features(
    adata=sds,
    n_features=embed_params["nfeat"],
    filter_lower_quantile=0.005,
    filter_upper_quantile=0.005,
    whitelist=None,
    blacklist=Path(blacklist_file),
    max_iter=1,
    inplace=True,
)

logger.info("Perform snapatac2 embedding.")
sa2.tl.spectral(
    adata=sds,
    n_comps=embed_params["ncomp"],
    features="selected",
    random_state=0,
    sample_size=None,
    distance_metric="cosine",
    weighted_by_sd=True,
    inplace=True,
)

# ** knn
logger.info(f"Use {knn_nm} for knn")
sa2.pp.knn(
    adata=sds,
    n_neighbors=knn_params["n"],
    use_dims=None,
    use_rep="X_spectral",
    method=knn_params["method"],
    inplace=True,
    random_state=0,
)


def run_leiden(r: float) -> np.ndarray:
    logger.info(f"Resolution: {r}")
    res = sa2.tl.leiden(
        adata=sds,
        resolution=r,
        objective_function=leiden_params["obj"],
        min_cluster_size= leiden_params["min_size"],
        random_state=leiden_params["seed"],
        use_leidenalg=False,
        weighted=leiden_params["weight"],
        inplace=False,
    )
    return res  # pyright: ignore


min_r: float = leiden_params["minr"]
max_r: float = leiden_params["maxr"]
step: float = leiden_params["byr"]
logger.info(f"Run resolutions from {min_r} to {max_r} with {step}.")
rs = np.arange(min_r, max_r + step, step)
rs = np.array([round(i, 1) for i in rs])
clusters = list(map(run_leiden, rs))
# n_cell by n_resolution
clustering_sum["leiden"] = np.stack(clusters, axis = 1)
clustering_sum["leiden_r"] = rs

# *** get silhouette score
def get_silhouetee(r: float) -> float:
    return cal_silhouette(
        adata=sds,
        resolution=r,
        use_rep="X_spectral",
        obj=leiden_params["obj"],
        rerun_leiden=True,
        n_sample=leiden_params["n_sample"],
        random_state=leiden_params["seed"],
        idx= int(np.where(rs == r)[0]),
        logger = logger
    )


if sds.X.shape[0] > leiden_params["min_size"]:
    logger.info("Calculate Silhouetee Score.")
    silhouettes = list(map(get_silhouetee, rs))
    clustering_sum["sihouettes"] = silhouettes

# *** other metrics: PAC and disp
nr = leiden_params["repeat"]
if nr > 1 and sds.X.shape[0] > leiden_params["min_size"]:
    logger.info(f"Get metrics of PAC and disp based on {nr} repeats.")

    def get_cemba(r: float) -> Tuple[float, float]:
        return cemba(
            adata=sds,
            resolution=r,
            objective_function=leiden_params["obj"],
            n_repeat=leiden_params["repeat"],
            n_sample=leiden_params["n_sample"],
        )

    metrics_cemba = list(map(get_cemba, rs))
    clustering_sum["cemba_disp"] = [a0 for a0, _ in metrics_cemba]
    clustering_sum["cemba_PAC"] = [a1 for _, a1 in metrics_cemba]


# ** run UMAP
logger.info("Get umap.")
umap(
    adata=sds,
    key_added="umap",
    random_state=umap_params["seed"],
    inplace=True,
    a=umap_params["a"],
    b=umap_params["b"],
    init=umap_params["init"],
    min_dist=umap_params["min_dist"],
    n_neighbors=umap_params["n_neigh"],
    metric=umap_params["metric"]
)

logger.info("Save embed and umap to clustering_sum object.")
clustering_sum['barcode'] = sds.obs_names
clustering_sum["embed"] = sds.obsm["X_spectral"]
clustering_sum["embed_eigs"] = sds.uns['spectral_eigenvalue']
clustering_sum["umap"] = sds.obsm["X_umap"]

logger.info(f"Dump clustering_sum object into {out_sum_file} using pickle.")
with open(out_sum_file, 'wb') as f:
    pickle.dump(clustering_sum, f)

# done
sds.close()
logger.info("Run leiden done.")
