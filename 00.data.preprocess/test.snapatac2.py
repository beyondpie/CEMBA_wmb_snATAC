import snapatac2 as snap
import sys
import pyprojroot
proj_root = pyprojroot.here()
packdir = f"{proj_root}/package/python"
sys.path.insert(0, packdir)
from utils import printfn
import h5py
from pathlib import Path
# import gc
# gc.collect()

import numpy as np
import argparse
import os

# parser = argparse.ArgumentParser("snapatac2 from bam")
# parser.add_argument("-b", "--bam", type = str,
#                     default = "bams/CEMBA190423_10A.bam")
# args = parser.parse_args()

# * load chrom size
with open(f"{proj_root}/meta/mm10.chrom.sizes.lite", 'r') as f:
    chrom_size = {l.strip().split()[0]: int(l.strip().split()[1])
                   for l in f.readlines()}
# gtf23 = f"{proj_root}/meta/gencode.vM23.primary.annot2.gtf"
gtf23 = f"{proj_root}/meta/gencode.vM23.gene.annot2.gtf"

bam_file = "./bams/CEMBA190423_10A.bam"
outdir = "fragments"
outfnm = "CEMBA190423_10A.frag"
frag_file = f"{outdir}/{outfnm}"

# * generate fragment_file from bam_file
# fragment has been corrected
stat = snap.pp.make_fragment_file(bam_file = Path(bam_file),
                                  output_file = Path(frag_file),
                                  is_paired = True,
                                  barcode_regex = "^(\w+):.+",
                                  shift_left = 4,
                                  shift_right = -5,
                                  min_mapq = 30)

## save stat to some places
# * get anndata from fragment
ann_file = f"{outdir}/CEMBA190423_10A.h5ad"

## this will inplace filter the cells
snap_data = snap.pp.import_data(fragment_file = Path(frag_file),
                                genome = snap.genome.mm10,
                                ## whitelist: file or list of barcode
                                whitelist = None,
                                ## output h5ad file
                                file = Path(ann_file),
                                min_num_fragments = 1,
                                min_tsse = 0,
                                shift_left = 0,
                                shift_right = 0,
                                chrom_size = chrom_size,
                                sorted_by_barcode = False)
snap_data.close()
snap_data = snap.read(filename = ann_file)
# snap_data.close() before rerun command above
# snap.pl.tsse(snap_data, interactive = True,
#              min_fragment = 1,
#              width = 500,
#              height = 500,
#              show = True,
#              out_file = None,
#              scale = None)

# * QC
index_qc = snap.pp.filter_cells(data = snap_data,
                                min_counts = 1000,
                                min_tsse = 10,
                                inplace = False)
snap_data.subset(obs_indices = index_qc,
                 out = Path("test.qc.hdf5"))
snap_data.close()

snap_data = snap.read(filename = ann_file)

# * add gene matrix
gmat_file = f"{outdir}/CEMBA190423_10A.gmat3.h5ad"
snap_gmat = snap.pp.make_gene_matrix(adata = snap_data,
                                     gene_anno = Path(gtf23),
                                     # gene_anno = snap.genome.mm10,
                                     file = gmat_file,
                                     id_type = "gene")

# * doublet removal
snap.pp.scrublet(snap_gmat, features = None)
snap_gmat.close()
snap_gmat = snap.read(filename = gmat_file)
snap_data2 = snap.pp.filter_doublets(adata = snap_gmat, inplace = True)

# * add bmat
snap.pp.add_tile_matrix(adata = snap_data,
                        bin_size = 500,
                        inplace = False,
                        exclude_chroms = ['chrM', 'Un', 'random'],
                        file = Path("test.bmat.hdf5"))

# * test gtf file
import pyranges as pr
a = pr.read_gtf(gtf23)
