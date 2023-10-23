import os
import sys
from pathlib import Path
import logging
import numpy as np
import numpy
import argparse
import shutil

import snapatac2 as sa2
import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils # type: ignore # noqa: E402

parser = argparse.ArgumentParser("snapatac2 L1 KNN")
parser.add_argument("--embed_file", type = str)
parser.add_argument("--outf", type = str, default = "test_knn.hdf5")
parser.add_argument("--kmethod", type = str, default = 'exact')
parser.add_argument("--knn", type = int, default = 50)
parser.add_argument("--logfile", type = str,
                    default = "log/test_sa2_l1_embed_knn.log")
parser.add_argument("--debug", type = int, default = 1)
parser.add_argument("-i", "--ipython", action = "store_true")
parser.add_argument("--simple-prompt", action = "store_true")

args = parser.parse_args()

# * set log
logger = utils.set_file_logger(fnm = args.logfile, #type: ignore
                               name = "sa2.04.l1_knn")
if args.debug == 0:
    debug = False
else:
    debug = True
    logger.warning("DEBUG mode is open")

# * meta
k = args.knn
km = args.kmethod
embed_file = args.embed_file
if not os.path.exists(embed_file):
    err_msg = f"{embed_file} is not found"
    logging.error(err_msg)
    sys.exit(err_msg)
outf = args.outf
outdir = os.path.dirname(outf)
if os.path.exists(outf):
    logger.warning(f"{outf} exists and remove it.")
    os.remove(outf)
else:
    os.makedirs(outdir, exist_ok = True)


# * main
logger.info(f"Copy {embed_file} to {outf}")
shutil.copyfile(src = embed_file,
                dst = outf)

logger.info(f"Read AnnData to RAM from : {outf}")
sds = sa2.read(outf, backend = None)

logger.info("Start to run KNN.")
sa2.pp.knn(
    adata = sds,
    n_neighbors = k,
    use_dims = None,
    use_rep = 'X_spectral',
    method = km,
    inplace = True,
    random_state = 0
)
sds.close()
logger.info("Done")
