import os
import sys
from pathlib import Path
import logging
import numpy as np
import tracemalloc
import argparse

import snapatac2 as sa2

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils # type: ignore # noqa: E402

parser = argparse.ArgumentParser("snapatac2 L1 embedding")
parser.add_argument("--samples", type = str,
                    default = "test_samplex.txt")
parser.add_argument("--bmatdir", type = str,
                    default = "snapatac2_pp_out/bmat")
parser.add_argument("--nfeature", type = int,
                    default = 500000)
parser.add_argument("--outf", type = str,
                    default = "sa2_l1_embed_500000_1.0_50_cosine.hdf5")
parser.add_argument("--blacklist", type = str,
                    default = "../meta/mm10.chrom.sizes.lite")
parser.add_argument("--samplesize", type = float,
                    default = 1.0)
parser.add_argument("--ncomps", type = int, default = 50)
parser.add_argument("--dmetric", type = str, default = 'cosine')
parser.add_argument("--knnmethod", type = str, default = 'exact')
parser.add_argument("--knn", type = int, default = 50)
parser.add_argument("--logfile", type = str,
                    default = "log/test_sa2_l1_embed_knn.log")
parser.add_argument("--debug", type = int, default = 1)
parser.add_argument("-i", "--ipython", action = "store_true")
parser.add_argument("--simple-prompt", action = "store_true")

args = parser.parse_args()

# * prepare outfile
nf = args.nfeature
ss = args.samplesize
nc = args.ncomps
dm = args.dmetric

# * set log
logger = utils.set_file_logger(fnm = args.logfile, #type: ignore
                               name = "sa2.03.l1_embed")

if args.debug == 0:
    debug = False
else:
    debug = True
    logger.warning("DEBUG mode is open.")

# * meta
blacklist_file = Path(args.blacklist)
if not os.path.exists(blacklist_file):
    err_msg = f"{blacklist_file} is not found."
    logging.error(err_msg)
    sys.exit(err_msg)
bmat_dir = Path(args.bmatdir)
if not os.path.exists(bmat_dir):
    err_msg = f"{bmat_dir} is not found."
    logging.error(err_msg)
    sys.exit(err_msg)

outf = args.outf
outdir = os.path.dirname(outf)
os.makedirs(outdir, exist_ok = True)
if os.path.exists(outf):
    logger.warning(f"{outf} exists and remove it.")
    os.remove(outf)
logger.info(f"SnapATAC2 L1 embedding: {outf}")


# * main
tracemalloc.start()
logger.info("Load samples")
with open(args.samples, 'r') as f:
    samples = [line.strip() for line in f.readlines() if len(line) > 1]
files = [(s, Path(f"{bmat_dir}/{s}.sa2.bmat.hdf5")) for s in samples]
logger.info("Load AnnDataSet to cover the samples.")
sds = sa2.AnnDataSet(adatas = files, filename = outf)

# logger.info("Update AnnDataSet barcodes' names.")
# FIXME: error at tscc: numpy.core._exceptions._UFuncNoLoopError:
# ufunc 'add' did not contain a loop with signature matching
# types (dtype('<U16'), dtype('<U1')) -> None
# sds.obs_names = np.array(sds.obs['sample']) + "." + np.array(sds.obs_names)

logger.info("Filtering features.")
    
sa2.pp.select_features(
    adata = sds,
    n_features = nf,
    filter_lower_quantile = 0.005,
    filter_upper_quantile = 0.005,
    whitelist = None,
    blacklist = blacklist_file,
    max_iter = 1,
    inplace = True
)

logger.info("Spectral embedding.")
s1 = tracemalloc.take_snapshot()
logger.debug("RAM before spectral: ")
for stat in s1.statistics('lineno')[:10]:
    logger.debug(stat)

sa2.tl.spectral(
    adata = sds,
    n_comps = nc,
    features = 'selected',
    random_state = 0,
    sample_size = ss,
    sample_method = 'random',
    distance_metric = dm,
    weighted_by_sd = True,
    inplace = True
)
logger.info("Spectral done.")
s2 = tracemalloc.take_snapshot()
logger.debug("RAM after spectral: ")
for stat in s2.statistics('lineno')[:10]:
    logger.debug(stat)

logger.debug("Diff RAM in embed:")
for stat in s2.compare_to(s1, 'lineno')[:10]:
    logger.debug(stat)

sds.close()
logger.info("Done")







