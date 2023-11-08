import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
import pyprojroot
packdir = f"{pyprojroot.here()}/package/python"
sys.path.insert(0, packdir)
import utils
import mycelloracle


parser = argparse.ArgumentParser(description = "run gimmemotifs")
parser.add_argument('--pdcbedpe', type = str)
parser.add_argument('--outdir', type = str, default = "tfscan")
parser.add_argument('--threshold', type = int, default = 10)
parser.add_argument('--refgenome', type = str, default = "mm10")

args = parser.parse_args()

# * main
# to celloracle supported file
print(f"TFscan for {args.pdcbedpe}")
pdc_bedpe = pd.read_table(args.pdcbedpe, sep = "\t", header = None)
pdc_bedpe.columns = ["pchr", "pstart", "pend", "dchr",
                     "dstart", "dend", "pair", "pearson", "n1", "n2"]

peak_id = pdc_bedpe['dchr']  + '_' + pdc_bedpe['dstart'].map(str) + '_'  + pdc_bedpe['dend'].map(str)
gene_short_name = pdc_bedpe['pair'].apply(lambda x: x.split('|')[0])
pdc = pd.concat([peak_id, gene_short_name], axis = 1)
pdc.columns = ["peak_id", 'gene_short_name']
mycelloracle.check_peak_format(peaks_df = pdc, ref_genome = args.refgenome)

# * tfcan
tfi = ma.TFinfo(peak_data_frame = pdc, ref_genome = args.refgenome)
tfi.scan(fpr = 0.02, motifs = None, verbose = True)
# save result
prefix = os.path.basename(args.pdcbedpe).replace(".bedpe", "")
tfi.to_hdf5(file_path=f"{args.outdir}/{prefix}.celloracle.tfinfo")

# * filter
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=args.threshold)
tfi.make_TFinfo_dataframe_and_dictionary(verbose = True)
df_tfi = tfi.to_dataframe()
# save result
df_tfi.to_parquet(f"{args.outdir}/{prefix}.baseGRN.df.parquet")
print("Done")







