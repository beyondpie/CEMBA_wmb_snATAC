import snapatac2 as sa2
import os
import sys
import argparse
import numpy as np
import logging
from pathlib import Path
from typing import Dict
import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils # type: ignore # noqa: E402

parser = argparse.ArgumentParser("snapatac2 preprocessing from bam file.")
parser.add_argument("--sample", type = str, default = "CEMBA190423_10A")
parser.add_argument("--outdir", type = str, default = "snapatac2_pp_out")
parser.add_argument("--bamfile",
                    type = str,
                    default = "snapatac2_pp_out/raw_fragment/CEMBA190423_10A.sa2.frag.hdf5")
parser.add_argument("--chromfile", type = str,
                    default = "../meta/mm10.chrom.sizes.lite")
parser.add_argument("--gtfile", type = str,
                    default = "../meta/gencode.vM23.gene.annot2.rm.semicolon.gtf")
parser.add_argument("--blacklist", type = str,
                    default = "../meta/mm10.blacklist.bed")
parser.add_argument("--logfile", type = str,
                    default = "log/CEMBA190423_10A_snapatac2_pp.log")
parser.add_argument("--debug", type = int, default = 1)
parser.add_argument("-i", "--ipython", action = "store_true")
parser.add_argument("--simple-prompt", action = "store_true")

args = parser.parse_args()

# * prepare outfile
out_dir = args.outdir
sample = args.sample
frag_file = f"{out_dir}/raw_fragment/{sample}.sa2.frag.hdf5"
frag_stat_file = f"{out_dir}/raw_stat/{sample}.sa2.frag.stat"
raw_file = f"{out_dir}/raw_fragment/{sample}.sa2.raw.hdf5"
qc_file = f"{out_dir}/qc/{sample}.sa2.qc.hdf5"
dlt_file = f"{out_dir}/doublet/{sample}.sa2.dlt.hdf5"
bmat_file = f"{out_dir}/bmat/{sample}.sa2.bmat.hdf5"
gmat_file = f"{out_dir}/gmat/{sample}.sa2.gmat.hdf5"

# * set log
logger = utils.set_file_logger(fnm = args.logfile,
                               name = "sa2.01.pp")
logger.info(sample)

if args.debug == 0:
    debug = False
else:
    debug = True
    logger.warning("DEBUG mode is open.")

# * prepare meta
with open(args.chromfile, 'r') as f:
    tmp_lines = [l.strip() for l in f.readlines()] # noqa
    chrom_sizes:Dict[str, int] = {
        l.split()[0]: int(l.split()[1]) for l in tmp_lines} # noqa
bam_file = args.bamfile
gtf_file = args.gtfile
if not os.path.exists(bam_file):
    err_msg = f"{bam_file} is not found."
    logging.error(err_msg)
    sys.exit(err_msg)
# if not os.path.exists(gtf_file):
#     err_msg = f"{gtf_file} is not found."
#     logging.error(err_msg)
#     sys.exit(err_msg)

# * main
logger.info("Get raw fragment from bam file.")
if debug and os.path.exists(frag_file):
    logger.info(f"{frag_file} exist in debug mode.")
else:
    if os.path.exists(frag_file):
        logger.warning(f"{frag_file} exists and remove it.")
        os.remove(frag_file)
    stat = sa2.pp.make_fragment_file(
        bam_file = Path(bam_file),
        output_file = Path(frag_file),
        is_paired = True,
        barcode_regex = "^(\w+):.+",
        shift_left = 4,
        shift_right = -5,
        min_mapq = 30
    )
    logger.info("Save raw fragment's stat.")
    with open(frag_stat_file, 'w') as f:
        f.writelines(str(stat))

logger.info("Load the fragment data for QC.")
if os.path.exists(raw_file):
    logger.warning(f"{raw_file} exists and remove it.")
    os.remove(raw_file)
sd_raw = sa2.pp.import_data(
    fragment_file = Path(frag_file),
    genome = sa2.genome.mm10, # type: ignore
    whitelist = None,
    file = Path(raw_file),
    min_num_fragments = 1,
    min_tsse = 0,
    shift_left = 0,
    shift_right = 0,
    chrom_size = chrom_sizes,
    sorted_by_barcode = False
)
logger.info("Perform QC")
index_qc = sa2.pp.filter_cells(
    data = sd_raw,
    min_counts = 1000,
    min_tsse = 10,
    inplace = False
)
logger.info("Save fragments after QC.")
if os.path.exists(qc_file):
    logger.warning(f"{qc_file} exists and remove it.")
    os.remove(qc_file)
sd_raw.subset(
    obs_indices = index_qc,
    out = Path(qc_file)
)
sd_raw.close()

# logger.info("Add gene matrix for doublet removal.")
sd_qc = sa2.read(Path(qc_file))

if os.path.exists(dlt_file):
    logger.warning(f"{dlt_file} exists and remove it.")
    os.remove(dlt_file)

## gmat: gene body + 2kb upstream
sd_dlt = sa2.pp.make_gene_matrix(
    adata = sd_qc,
    gene_anno = sa2.genome.mm10, # type: ignore
    # gene_anno = Path(gtf_file),
    file = Path(dlt_file),
    use_x = False,
    id_type = "gene"
)
logger.info("Perform doublet removal based on gene matrix")
sa2.pp.scrublet(sd_dlt, features = None,
                n_comps = 15,
                sim_doublet_ratio = 2,
                expected_doublet_rate = 0.08,
                n_neighbors = None,
                use_approx_neighbors = True,
                random_state = 0,
                inplace = True)
sd_dlt.uns['reference_sequences'] = sd_qc.uns['reference_sequences']
sd_dlt.obsm['insertion'] = sd_qc.obsm['insertion']

logger.info("Filter doublets.")
sa2.pp.filter_doublets(
    adata = sd_dlt,
    probability_threshold = 0.5,
    score_threshold=None,
    inplace = True,
    verbose = True
)
sd_dlt.close()

logger.info("Get bmat for the samples and save it.")
sd_dlt = sa2.read(Path(dlt_file))
if os.path.exists(bmat_file):
    logger.warning(f"{bmat_file} exists and remove it.")
    os.remove(bmat_file)
sd_bmat = sa2.pp.add_tile_matrix(
    adata = sd_dlt,
    bin_size = 500,
    inplace = False,
    exclude_chroms= ['chrM', 'M', 'chrUn'],
    file = Path(bmat_file)
)
sd_bmat.close()
# logger.info("Get modified vM23 gmat for the samples.")
# sd_gmat = sa2.pp.make_gene_matrix(
#     adata = sd_dlt,
#     gene_anno = Path(gtf_file),
#     file = None,
#     id_type = "gene"
# )
# logger.info("Save modified vM23 gmat to file.")
# sd_gmat.write(filename = Path(gmat_file))
logger.info("Done.")
