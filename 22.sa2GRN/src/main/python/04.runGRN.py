import argparse
import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc
import pyarrow.parquet as pq
import pyprojroot

proj_dir = pyprojroot.here()
packdir = f"{proj_dir}/package/python"
sys.path.insert(0, packdir)
import utils
import mycelloracle

work_dir = os.path.join(proj_dir, "22.sa2GRN")
rsc_dir = os.path.join(work_dir, "src/main/resource")

parser = argparse.ArgumentParser(description = "run GRN")
parser.add_argument("--tfdir", type = str,
                    default = os.path.join(work_dir, "out/tfscan"))
parser.add_argument("--scRNAdir", type = str,
                    default = os.path.join(rsc_dir,
                                           "sa2.allen.logCPM.vf3281.ds1000.subclass.specific"))
parser.add_argument("--outdir", type = str,
                    default = os.path.join(work_dir, "out/GRN"))
parser.add_argument("--group", type = str, default = "TU-ARH_Otp_Six6_Gaba")
parser.add_argument("--groupcol", type = str, default = "subclass")
parser.add_argument("--baseGRN_all", type = float, default = 0)
parser.add_argument("--njobs", type = int, default = 2)
parser.add_argument("--pvalue", type = float, default = 0.001)
parser.add_argument("--nlink", type = int, default = 10000)
parser.add_argument("--alpha", type = int, default = 10)



args = parser.parse_args()
group = args.group
group_col = args.groupcol
njobs = args.njobs
alpha = args.alpha
kmin4impute = 4
## filter GRN
p = args.pvalue
ntoplink = args.nlink
outdir = args.outdir
os.makedirs(outdir, exist_ok = True)
baseGRN_all = args.baseGRN_all
# * load data
print(f"Loading data for {group}.")
adata = sc.read_h5ad(
    f"{args.scRNAdir}/sa2.allen.subclass.{group}.ann.hdf5")
if baseGRN_all > 0:
    print("Use global baseGRN from sa2subclass.all.baseGRN.df.parquet.")
    tfi_all = mycelloracle.load_tfidf(
        f"{args.tfdir}/sa2subclass.all.baseGRN.df.parquet")
else:
    print(f"Use {group} baseGRN from sa2subclass.all.baseGRN.df.parquet.")
    tfi_all = mycelloracle.load_tfidf(
        f"{args.tfdir}/sa2subclass.{group}.pdc.baseGRN.df.parquet"
    )
    
# mycelloracle.check_peak_format(tfi_all, ref_genome = "mm10")

# * scRNAseq preprocessing
print("Preprocessing scRNAseq.")
# sc.pp.filter_genes(adata, min_counts = 1)
# sc.pp.normalize_per_cell(adata, key_n_counts = 'n_counts_all')
# adata.raw = adata
## if no raw_count, oracle will raise error
adata.layers["raw_count"] = adata.X.copy()
# sc.pp.log1p(adata)
adata.layers["logCPM"] = adata.X.copy()
sc.pp.scale(adata)

# * pca
sc.tl.pca(adata, svd_solver = "arpack")
# used for umap
sc.pp.neighbors(adata, n_neighbors = 4, n_pcs = 50)
sc.tl.umap(adata)

# * create cell oracle
print("Creat CellOracle object.")
oracle = co.Oracle()
adata.X = adata.layers["logCPM"].copy()
oracle.import_anndata_as_normalized_count(adata = adata,
                                          cluster_column_name = group_col,
                                          embedding_name = 'X_umap')
oracle.import_TF_data(TF_info_matrix = tfi_all)
oracle.perform_PCA()
# plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(
    np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]
# plt.axvline(n_comps, c="k")
# plt.show()
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
k = max(int(0.025 * n_cell), kmin4impute)

print("KNN imputation.")
oracle.knn_imputation(n_pca_dims = n_comps, k = k, balanced = True,
                      b_sight = k * 8, b_maxl = k*4, n_jobs = njobs)
# get GRN
print(f"Infer the GRN for {group}")
links = oracle.get_links(cluster_name_for_GRN_unit = group_col,
                         alpha = alpha, verbose_level = 0,
                         n_jobs = njobs)
# use default parameter p and threshold_number
# p=0.001, weight='coef_abs', threshold_number=10000
links.filter_links(p=p, weight="coef_abs",
                   threshold_number=ntoplink)
# network analysis
links.get_network_score()

# visualize results
# links.plot_degree_distributions(plot_model=True)
# links.get_network_score()
# links.plot_scores_as_rank('MB-HB_Lhx1_Glut')
# links.merged_score.head()

# save GRN
print("Save results.")
links.to_hdf5(file_path=f"{outdir}/GRN.{group}.celloracle.links")
links.links_dict[group].to_csv(f"{outdir}/rawGRN.{group}.csv")
print("Done")
