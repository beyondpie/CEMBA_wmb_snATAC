import os
import sys
import traceback
from typing import Dict, Tuple
from pathlib import Path
import pickle

import numpy as np

from numba.core.errors import NumbaDeprecationWarning 
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import snapatac2 as sa2 # noqa: E402
import pyprojroot # noqa: E402

code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils  # noqa: E402
from leiden import cemba  # pyright: ignore # noqa: E401, E402, F401
from leiden import cal_silhouette  # pyright: ignore # noqa: E402, F401
from leiden import umap  # pyright: ignore # noqa: E402

# * log
logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "leiden.clustering"
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
leiden_params: Dict = snakemake.params["leiden"] #pyright: ignore # noqa
leiden_nm: str = leiden_params["name"]
snap_file = snakemake.input["snap_file"][0] #pyright: ignore # noqa
embed_params: Dict = snakemake.params["embed"]  # pyright: ignore # noqa
knn_params: Dict = snakemake.params["knn"]  # pyright: ignore # noqa
umap_params: Dict = snakemake.params["umap"] #pyright: ignore # noqa

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
logger.info(f"Run leiden with {leiden_nm}.")
logger.info("Load snapatac2 file.")
logger.info(f"Load snapatac2 anndataset from {snap_file}.")
sds = sa2.read(Path(snap_file))

# * run leiden with single/mulitple resolutions
def run_leiden(r:float) -> np.ndarray:
    logger.info(f"Resolution: {r}")
    res = sa2.tl.leiden(
        adata = sds,
        resolution = r,
        objective_function=leiden_params["obj"],
        min_cluster_size= leiden_params["min_size"],
        random_state = leiden_params["seed"],
        use_leidenalg = False,
        weighted = leiden_params["weight"],
        inplace = False
    )
    return res #pyright: ignore

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
            adata = sds,
            resolution = r,
            objective_function = leiden_params["obj"],
            n_repeat = leiden_params["repeat"],
            n_sample = leiden_params["n_sample"]
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


    
