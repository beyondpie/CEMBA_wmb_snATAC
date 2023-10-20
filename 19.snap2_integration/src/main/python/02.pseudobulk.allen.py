import os
import sys

import pandas as pd
import anndata as ad

from pyprojroot import here

proj_root = str(here())
sys.path.insert(0, os.path.join(proj_root, "package/python"))
from myanndata import grouped_obs_mean, get_barcode2group
from utils import set_file_logger

# * configs
resource_dir = f"{proj_root}/19.snap2_integration/src/main/resource"
logger = set_file_logger(fnm=f"{resource_dir}/avgexp.log", name="avgexp")
annot_allen = pd.read_table(
    f"{resource_dir}/AIT21_annotation.tsv", sep="\t", header=0, index_col=False
)

# get barcode annotation info for allen's data.
groups_allen = [
    "cluster_id",
    "supertype_id",
    "subclass_id",
    "class_id",
    "nt_type_label",
    "nt_type_combo_label",
]


# * handle allen neuron data
def get_avgexp_from_allen(
    ann: ad.AnnData,
    group_meta: pd.DataFrame,
    group_key: str = "subclass_id",
    out_dir: str | None = None,
    cell_class: str = "neuron",
    tech: str = "10xv3",
) -> pd.DataFrame:
    avgexp = grouped_obs_mean(adata=ann, group_meta=group_meta, group_key=group_key)
    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        avgexp.to_csv(
            os.path.join(out_dir, f"{cell_class}_{tech}_{group_key}.avgexp.csv")
        )
    return avgexp


# * main
for tech in ["10xv3", "10xv2", "multiome"]:
    for cell_class in ["nn", "neuron"]:
        fnm = f"{resource_dir}/{cell_class}_male_{tech}_ann.h5ad"
        if not os.path.exists(fnm):
            print(f"{fnm} does not found.")
            continue
        logger.info(f"Loading data from: {fnm}")
        ann_data = ad.read_h5ad(filename=fnm, backed=None)
        logger.info("Finish loading file.")
        b2g = get_barcode2group(
            obs=ann_data.obs, annot=annot_allen, fromcol="cl", tocol=groups_allen
        )
        logger.info("start to get avgexp.")
        for g in groups_allen:
            tmp = get_avgexp_from_allen(
                ann_data,
                b2g,
                g,
                f"{resource_dir}/rknn",
                cell_class=cell_class,
                tech=tech,
            )


# ** merge 10xv3, 10xv2, multiome
def load_data(cell_class: str = "nn", tech="10xv3") -> ad.AnnData:
    fnm = os.path.join(resource_dir, f"{cell_class}_male_{tech}_ann.h5ad")
    return ad.read_h5ad(filename=fnm, backed=None)


nn_10xv3 = load_data("nn", "10xv3")
neuron_10xv3 = load_data("neuron", "10xv3")
nn_10xv2 = load_data("nn", "10xv2")
neuron_10xv2 = load_data("neuron", "10xv2")
neuron_multiome = load_data("neuron", "multiome")

# validate obs_names are unique
# generate codes for get obs_names

# check obs_name is unique along the datasets
# 2263577 and they are unique
# all_obs_names = nn_10xv3.obs_names.append(
#     neuron_10xv3.obs_names).append(
#         nn_10xv2.obs_names).append(
#             neuron_10xv2.obs_names).append(
#                 neuron_multiome.obs_names)

# merge them
# cost near 800G RAM during copy
# after that about 500G RAM in memory
# So all_merged should be about 250G RAM
all_merged = ad.concat(
    [nn_10xv3, neuron_10xv3, nn_10xv2, neuron_10xv2, neuron_multiome]
)

b2g = get_barcode2group(
    obs=all_merged.obs, annot=annot_allen, fromcol="cl", tocol=groups_allen
)
b2g["barcode"] = b2g.index
b2g.to_csv(f"{resource_dir}/allen_male_all.b2g.csv", index=False)

# during this process, it cost about 700G RAM at most
for g in groups_allen:
    print(f"current group: {g}")
    tmp = get_avgexp_from_allen(
        all_merged,
        b2g,
        g,
        f"{resource_dir}/rknn",
        cell_class="all",
        tech="10x_multiome",
    )
# save all_merged
# check if some columns has different types
# Function to get the type of each element
def find_type(val):
    return type(val)

# Apply the function to the DataFrame
type_df = all_merged.obs.applymap(find_type)
# Check if any column has more than one unique type
mixed_columns = type_df.nunique() > 1
# Print columns with mixed types
print(mixed_columns[mixed_columns])
# delete the mixed type columns
mixed_columns_names = mixed_columns.index[mixed_columns]
## an external source column is mixed type
df_cleaned = all_merged.obs.drop(columns=mixed_columns_names)
all_merged.obs = df_cleaned

all_merged.write(f"{resource_dir}/allen_male_all.h5ad")
