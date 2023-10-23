import os
import sys
import importlib
import pandas as pd
import anndata as ad

from pyprojroot import here
proj_root: str = str(here())
sys.path.insert(0, os.path.join(proj_root, "package/python"))
from utils import set_file_logger
import cembav2env
importlib.reload(cembav2env)

# * configs
rsc_dir = os.path.join(proj_root,
                       "19.snap2_integration/"
                       "src/main/resource")
allen = cembav2env.Allen()
sa2atac = cembav2env.Sa2ATAC()

# * reduce anndata by removing raw data
# run only once
# let's keep using log-norm data for later integration
# since scaling is invariant towards normalization constant.

reduced_allen_dir = os.path.join(rsc_dir, "norawdata_allen")
if not os.path.exists(reduced_allen_dir):
    os.makedirs(reduced_allen_dir, exist_ok = False)
ann_10xv3_nn = allen.get_10xv3_nn_ann()
allen.get_reduced_ann(
    ann_10xv3_nn,
    outfnm = os.path.join(reduced_allen_dir, "nn_male_10xv3_ann_noraw.h5ad"))
del ann_10xv3_nn

ann_10xv3_neuron = allen.get_10xv3_neuron_ann()
allen.get_reduced_ann(
    ann_10xv3_neuron,
    outfnm = os.path.join(reduced_allen_dir,
                          "neuron_male_10xv3_ann_noraw.h5ad"))
del ann_10xv3_neuron


ann_all_nn = allen.get_allen_nn_ann()
allen.get_reduced_ann(
    ann_all_nn,
    outfnm = os.path.join(reduced_allen_dir,
                          "nn_male_allen_ann_noraw.h5ad"))
del ann_all_nn

ann_all_neuron = allen.get_allen_neuron_ann()
allen.get_reduced_ann(
    ann_all_neuron,
    outfnm = os.path.join(reduced_allen_dir,
                          "neuron_male_allen_ann_noraw.h5ad"))
del ann_all_neuron

# * remove raw from our atac gmat
# run only once
reduced_atac_dir = os.path.join(rsc_dir, "norawdata_atac")
if not os.path.exists(reduced_atac_dir):
    os.makedirs(reduced_atac_dir, exist_ok=True)
atac_ann = sa2atac.load_sa2gmat_ann()
barcode2L3: pd.DataFrame = sa2atac.read_barcode2L3()
sa2atac.add_L3_to_atac_ann(atac_ann, barcode2L3)


atac_rough_annot = sa2atac.read_rough_annot()
barcode2annot = atac_rough_annot.loc[barcode2L3["L3"]]
barcode2annot.index = barcode2L3["barcode"]

atac_nn_ann = sa2atac.get_nn_atac_ann(
    atac_ann, barcode2annot,
    outfnm = os.path.join(reduced_atac_dir, "nn_gmat_atac_ann.h5ad"))
del atac_nn_ann

atac_neuron_ann = sa2atac.get_neuron_atac_ann(
    atac_ann, barcode2annot,
    outfnm = os.path.join(reduced_atac_dir, "neuron_gmat_atac_ann.h5ad")
)
del atac_neuron_ann








