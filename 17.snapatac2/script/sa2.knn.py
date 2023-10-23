import os
from pathlib import Path
import sys
import traceback
from typing import Dict

import snapatac2 as sa2
import pyprojroot

code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils #pyright: ignore # noqa: E402

# * log
logger = utils.set_file_logger( #pyright: ignore # noqa
    fnm = snakemake.log[0], #pyright:ignore # noqa: F821
    name = "sa2.L2.embed"
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
snap_file = snakemake.input["snap_file"][0] #pyright: ignore # noqa
knn_params: Dict = snakemake.params["knn"] #pyright: ignore # noqa 
knn_nm: str = knn_params["name"]
logger.info(f"Use {knn_nm} for knn.")
logger.info(f"Load snapatac2 anndataset: {snap_file}.")
sds = sa2.read(Path(snap_file))
sa2.pp.knn(
    adata = sds,
    n_neighbors = knn_params["n"],
    use_dims = None,
    use_rep = "X_spectral",
    method = knn_params["method"],
    inplace = True,
    random_state = 0
)
sds.close()
logger.info("Run KNN done.")



