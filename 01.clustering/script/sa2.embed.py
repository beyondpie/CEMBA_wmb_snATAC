import os
import sys
import traceback
from pathlib import Path
from typing import Dict
import shutil
import pandas as pd
import numpy as np

# import snakemake
import snapatac2 as sa2
import pyprojroot

code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils #pyright: ignore # noqa: E402

# * log
logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "sa2.embed"
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
snap_file = snakemake.output["snap_file"][0]  #pyright: ignore # noqa
embed_params: Dict = snakemake.params['embed'] #pyright: ignore # noqa
embed_nm:str = embed_params["name"]

blacklist_file = snakemake.input['blacklist'] #pyright: ignore # noqa
barcode2id_file = snakemake.input["barcode2id"]  # pyright: ignore # noqa
cemba_anndata_file = snakemake.input["cemba_anndata"]  # pyright: ignore # noqa

# * main
# ** load cemba anndataset
logger.info(f"Loading cemba all AnnData: {cemba_anndata_file}.")
# use read only mode for reading the same data simutaneously
# sds_all = sa2.read_dataset(Path(cemba_anndataset_file), mode = 'r')
sds_all = sa2.read(Path(cemba_anndata_file), 'r')

# clustering info
cid = snakemake.params["cluster_id"]  # pyright: ignore # noqa
clevel = snakemake.params["clevel"]  # pyright: ignore # noqa

# ** subset sds_all
logger.info(f"Subset above AnnDataSet: {clevel}_{cid}.")
barcode2id: pd.DataFrame = pd.read_csv(barcode2id_file,
                                       header=0, index_col=False,
                                       dtype = str)
# cid here is treated as str
sub_barcodes = barcode2id[barcode2id[clevel] == cid]['barcode']
logger.info(f"Under {clevel}_{cid}, find {sub_barcodes.size} barcodes.")

# sds = sds_all.to_adata(
#     obs_indices = is_in_sub,
#     copy_x = True,
#     file = snap_file
# )

if sds_all.shape[0] == sub_barcodes.size:
    logger.info("Subset cells have the same size, so will only copy file.")
    sds_all.close()
    shutil.copyfile(cemba_anndata_file, snap_file)
    sds = sa2.read(snap_file, 'r+')
else:
    logger.info(f"Use to_data to subet the sds_all to {snap_file}.")
    all_barcodes = sds_all.obs_names
    # cost lots of time when sub_barcodes and all_barcodes are large.
    # is_in_sub = np.isin(all_barcodes, sub_barcodes, assume_unique = True)
    # NOTE: alternate np.isin
    # https://github.com/numpy/numpy/issues/14997
    a = set(sub_barcodes)
    is_in_sub = np.array([b in a for b in all_barcodes])
    logger.info(f"Find {is_in_sub.sum()} barcodes in cemba.")
    sds = sds_all.subset(
        obs_indices = is_in_sub,
        out = snap_file)

logger.info("Select features.")
sa2.pp.select_features(
    adata = sds,
    n_features = embed_params["nfeat"],
    filter_lower_quantile = 0.005,
    filter_upper_quantile = 0.005,
    whitelist = None,
    blacklist = blacklist_file,
    max_iter = 1,
    inplace = True
)

logger.info("Perform snapatac2 embedding.")
sa2.tl.spectral(
    adata = sds,
    n_comps = embed_params["ncomp"],
    features = "selected",
    random_state = 0,
    sample_size = embed_params["nsample"],
    distance_metric = "cosine",
    weighted_by_sd = True,
    inplace = True
)

sds.close()
logger.info("Run embed done")




