import snapatac2 as sa2
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils # type: ignore # noqa: E402
from utils import head, printfn # type: ignore # noqa: E402


# * load embed and knn
# ** all the features, exact knn
sds_afek = sa2.read("snapatac2_pp_out/l1ek/l1_ek_1.0_50_cosine_exact_50.hdf5")
# sds obs_names does not contain sample, then add it.
sds_afek.obs_names = np.array(sds_afek.obs['sample']) \
    + '.' \
    + np.array(sds_afek.obs_names)

# clustering
# about 10-20 minutes
sa2.tl.leiden(
    adata = sds_afek,
    resolution = 1.0,
    objective_function = 'modularity',
    random_state = 0,
    key_added = 'leiden',
    use_leidenalg = False,
    weighted = False,
    inplace = True
)

# umap
# about 30 minutes
sa2.tl.umap(adata = sds_afek,
            n_comps = 2,
            use_rep = 'X_spectral',
            key_added = 'umap',
            inplace = True)

# visualize
outf_umap = "L1_UMAP_allfea_50nc_ekm.pdf"
umap_plt = sa2.pl.umap(adata = sds_afek,
            color = 'leiden',
            sample_size = None,
            interactive = False,
            show = True,
            out_file = outf_umap)
sds_afek.close()

# ** top 500,000 features
sds = sa2.read(os.path.join("snapatac2_pp_out/l1_silencer",
                            "sa2_l1knn_500000_50_1.0_cosine_exact_50.hdf5"))
sds.obs_names = np.array(sds.obs['sample']) \
    + '.' \
    + np.array(sds.obs_names)
sa2.tl.leiden(
    adata = sds,
    resolution = 1.0,
    objective_function = 'modularity',
    random_state = 0,
    key_added = 'leiden',
    use_leidenalg = False,
    weighted = False,
    inplace = True
)
sa2.tl.umap(adata = sds,
            n_comps = 2,
            use_rep = 'X_spectral',
            key_added = 'umap',
            inplace = True)
outf_umap = os.path.join("snapatac2_pp_out",
                         "L1_UMAP_top0.5mfea_50nc_ekm.pdf")
umap_plt = sa2.pl.umap(adata = sds,
            color = 'leiden',
            sample_size = None,
            interactive = False,
            show = True,
            out_file = outf_umap)
sds.close()

# * load two modalities, check consensus
sds_all = sa2.read(os.path.join(
    "snapatac2_pp_out/l1ek",
    "l1_ek_1.0_50_cosine_exact_50.hdf5"
))
barcodes_sdsall = sds_all.obs_names
leiden_sdsall = sds_all.obs['leiden']
umap_sdsall = pl.concat(
    [pl.DataFrame({"barcodes": barcodes_sdsall}),
     pl.from_numpy(sds_all.obsm['X_umap'],
                   schema = ["UMAP1", "UMAP2"])],
    how = 'horizontal'
)

umap_sdsall.write_csv(
    file = os.path.join("snapatac2_pp_out/eval_l1",
                        "L1_UMAP_allfeat.csv"),
    has_header = True
)

L1_sdsall = pl.DataFrame(
    {"barcodes": barcodes_sdsall,
     "leiden_all" : leiden_sdsall })
sds_all.close()

sds_top = sa2.read(os.path.join(
    "snapatac2_pp_out/l1_silencer",
    "sa2_l1knn_500000_50_1.0_cosine_exact_50.hdf5"
))
barcodes_sdstop = sds_top.obs_names
leiden_sdstop = sds_top.obs['leiden']

umap_sdstop = pl.concat(
    [pl.DataFrame({"barcodes": barcodes_sdstop}),
     pl.from_numpy(sds_top.obsm['X_umap'],
                   schema = ["UMAP1", "UMAP2"])],
    how = 'horizontal'
)
umap_sdstop.write_csv(
    file = os.path.join("snapatac2_pp_out/eval_l1",
                        "L1_UMAP_topfeat.csv"),
    has_header = True
)


L1_sdstop = pl.DataFrame(
    {"barcodes": barcodes_sdstop,
     "leiden_top": leiden_sdstop
    }
)
sds_top.close()

L1_sds =  L1_sdsall.join(L1_sdstop, on = "barcodes")
# save results for annotations
L1_sds.write_csv(
    file = os.path.join("snapatac2_pp_out/eval_l1",
                        "L1_sds.csv"),
    has_header = True
)










