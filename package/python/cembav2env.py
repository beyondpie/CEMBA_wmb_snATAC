import os
from typing import List, Dict
import sys
import pandas as pd
import numpy as np
import anndata as ad
from pyprojroot import here

proj_dir: str = str(here())
python_dir: str = f"{proj_dir}/package/python"

class Sa2ATAC(object):
    def __init__(self) -> None:
        self.proj_root: str = str(here())
        self.intgn_rsc_dir: str = os.path.join(
            self.proj_root, "19.snap2_integration",
            "src/main/resource", "sa2_gmat"
        )
        self.noraw_dir = os.path.join(
            self.proj_root, "19.snap2_integration",
            "src/main/resource","norawdata_atac")
        self.nn_gmat_lognorm_ann = os.path.join(
            self.noraw_dir, "nn_gmat_atac_ann.h5ad")
        self.neuron_gmat_lognorm_ann = os.path.join(
            self.noraw_dir, "neuron_gmat_atac_ann.h5ad")
        self.sa2gmat_fnm: str = os.path.join(
            self.intgn_rsc_dir,
            "sa2_sa2default_ann_with_CPMlognorm.h5ad"
        )
        self.rough_annot_fnm: str = os.path.join(
            self.proj_root,
            "19.snap2_integration",
            "src/main/resource",
            "sa2_dlt2_L3toSa1Annot.rough.csv"
        )
        self.barcode2L3_fnm = os.path.join(
            self.proj_root, "17.snapatac2", "post_script",
            "sa2.barcode2L3.csv")
        self.barcode2clu2annot = os.path.join(
            self.proj_root, "19.snap2_integration",
            "src/main/resource",
            "sa2_dlt2_L3toSa1Annot.rough.csv")

    def read_nn_gmat_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.nn_gmat_lognorm_ann)
    def read_neuron_gmat_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_gmat_lognorm_ann)
    
    def load_sa2gmat_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.sa2gmat_fnm)
    def read_rough_annot(self) -> pd.DataFrame:
        r = pd.read_csv(self.rough_annot_fnm, header=0, sep=",")
        r.set_index("clu.sa2", drop=False, inplace = True)
        return r
    
    def read_barcode2L3(self) -> pd.DataFrame:
        r = pd.read_csv(self.barcode2L3_fnm, sep = ",", header = 0)
        r.set_index("barcode", drop = False, inplace = True)
        return r

    @staticmethod
    def subset_atac_ann(
        ann: ad.AnnData,
        meta: pd.DataFrame,
        which_col: str = "main_class",
        keys: List[str] = ["NN"],
    ) -> ad.AnnData:
        if which_col not in meta.columns:
            raise KeyError(f"{which_col} does not exist.")
        if not meta[which_col].isin(keys).any():
            raise KeyError("keys do not exist in column {which_col}.")
        if not meta.index.equals(ann.obs.index):
            meta = meta.loc[ann.obs_names]
        return ann[meta[which_col].isin(keys)]

    @staticmethod
    def generate_nn_atac_ann(ann: ad.AnnData, rough_annot: pd.DataFrame,
                        outfnm: str|None = None) -> ad.AnnData:
        r = Sa2ATAC.subset_atac_ann(ann, rough_annot, "mainclass", ["NN"])
        r.layers = None
        if outfnm is not None:
            r.write(outfnm)
        return r 

    @staticmethod
    def generate_neuron_atac_ann(ann: ad.AnnData,
                            rough_annot: pd.DataFrame,
                            outfnm:str|None = None) -> ad.AnnData:
        r = Sa2ATAC.subset_atac_ann(ann,
                                    rough_annot,
                                    "mainclass", ["GLUT", "GABA"])
        r.layers = None
        if outfnm is not None:
            r.write(outfnm)
        return r
    
    @staticmethod
    def add_L3_to_atac_ann(ann: ad.AnnData,
                           barcode2L3: pd.DataFrame) -> None:
        barcode2L3.set_index("barcode", drop = False)
        if not ann.obs.index.equals(barcode2L3.index):
            barcode2L3.loc[ann.obs.index]
        ann.obs["L3"] = barcode2L3["L3"]
        return None

class Allen(object):
    def __init__(self) -> None:
        self.proj_root: str = str(here())
        self.intgn_rsc_dir: str = os.path.join(
            self.proj_root, "19.snap2_integration", "src/main/resource"
        )
        self.nn_10xv3_ann_fnm: str = os.path.join(
            self.intgn_rsc_dir, "nn_male_10xv3_ann.h5ad"
        )
        self.nn_10xv2_ann_fnm: str = os.path.join(
            self.intgn_rsc_dir, "nn_male_10xv2_ann.h5ad"
        )
        self.neuron_10xv3_ann_fnm: str = os.path.join(
            self.intgn_rsc_dir, "neuron_male_10xv3_ann.h5ad"
        )
        self.neuron_10xv2_ann_fnm: str = os.path.join(
            self.intgn_rsc_dir, "neuron_male_10xv2_ann.h5ad"
        )
        self.neuron_male_sn_fnm: str = os.path.join(
            self.intgn_rsc_dir, "neuron_male_multiome_ann.h5ad"
        )
        self.allen_male_all_ann_fnm: str = os.path.join(
            self.intgn_rsc_dir, "allen_male_all.h5ad"
        )
        self.noraw_dir: str = os.path.join(self.intgn_rsc_dir,
                                           "norawdata_allen")
        self.nn_10xv3_lognorm_fnm = os.path.join(self.noraw_dir,
                                         "nn_male_10xv3_ann_noraw.h5ad")
        self.nn_allen_lognorm_fnm = os.path.join(self.noraw_dir,
                                       "nn_male_allen_ann_noraw.h5ad")
        self.neuron_10xv3_lognorm_fnm = os.path.join(
            self.noraw_dir,"neuron_male_10xv3_ann_noraw.h5ad")
        self.neuron_allen_lognorm_fnm = os.path.join(
            self.noraw_dir, "neuron_male_allen_ann_noraw.h5ad")

    def read_nn_10xv3_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.nn_10xv3_lognorm_fnm)
    def read_nn_allen_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.nn_allen_lognorm_fnm)
    def read_neuron_10xv3_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_10xv3_lognorm_fnm)
    def read_neuron_allen_lognorm_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_allen_lognorm_fnm)

    # all the later get_[10xv3|10xv2|all|sn]_[nn|neuron]_ann function
    # load the original ann, having raw and log-normalized data
    def get_10xv3_nn_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.nn_10xv3_ann_fnm)

    def get_10xv2_nn_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.nn_10xv2_ann_fnm)

    def get_10xv3_neuron_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_10xv3_ann_fnm)

    def get_10xv2_neuron_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_10xv2_ann_fnm)

    def get_sn_neuron_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.neuron_male_sn_fnm)

    def get_allen_all_ann(self) -> ad.AnnData:
        return ad.read_h5ad(self.allen_male_all_ann_fnm)

    def get_allen_neuron_ann(self) -> ad.AnnData:
        neuron_10xv3 = self.get_10xv3_neuron_ann()
        neuron_10xv2 = self.get_10xv2_neuron_ann()
        neuron_sn = self.get_sn_neuron_ann()
        return ad.concat([neuron_10xv3, neuron_10xv2, neuron_sn])

    def get_allen_nn_ann(self) -> ad.AnnData:
        nn_10xv3 = self.get_10xv3_nn_ann()
        nn_10xv2 = self.get_10xv2_nn_ann()
        return ad.concat([nn_10xv3, nn_10xv2])

    @staticmethod
    def grouped_obs_mean(
        adata: ad.AnnData,
        group_meta: pd.DataFrame,
        group_key: str,
        layer: str | None = None,
        gene_symbols: List[str] | None = None,
    ) -> pd.DataFrame:
        """
        Ref: https://github.com/scverse/scanpy/issues/181
        Add checks before running the content.
        """
        # check
        if (layer is not None) and (layer not in adata.layers.keys()):
            raise KeyError(f"{layer} not found in adata layer.")
        # group_meta should have two columns: barcode, group_key
        if ("barcode" not in group_meta.columns) or (
            group_key not in group_meta.columns
        ):
            raise KeyError(f"Either barcode or {group_key} not in group_meta.")
        # and group_meta index by barcode
        if adata.shape[0] != group_meta.shape[0]:
            raise RuntimeError("adata and group_meta have diff # of cells.")

        # make sure the anndata and group_meta have the same order of barcodes
        group_meta.set_index("barcode", drop=False)
        group_meta = group_meta.loc[adata.obs_names]

        idy = [True] if gene_symbols is None else adata.var_names.isin(gene_symbols)
        if sum(idy) == 0:
            raise KeyError("No gene_symbols found in adata.")
        match gene_symbols:
            case None:
                subset_adata = adata
            case _:
                print(f"{sum(idy)} of {len(gene_symbols)}",
                      "genes found in adata.")
                subset_adata = adata[:, idy]
        match layer:
            case None:
                X = subset_adata.X
            case _:
                print(f"use {layer} for X.")
                X = subset_adata.layers[layer]
        grouped = group_meta.groupby(group_key)
        print(f"find {grouped.ngroups} under group {group_key}")
        out = pd.DataFrame(
            np.zeros((subset_adata.shape[1], len(grouped)), dtype=np.float64),
            columns=list(grouped.groups.keys()),
            index=subset_adata.var_names,
        )
        for group, idx in grouped.indices.items():
            out[group] = np.ravel(X[idx, :].mean(axis=0, dtype=np.float64))
        return out

    @staticmethod
    def check_pandas_column_type_consistent(obs: pd.DataFrame) -> pd.DataFrame:
        type_df = obs.applymap(lambda x: type(x))
        mixed_columns = type_df.nunique() > 1
        if sum(mixed_columns) > 0:
            mixed_column_names = mixed_columns.index[mixed_columns]
            print(f"find inconsistent column: {mixed_column_names}.")
            print("will remove from the obs.")
            return obs.drop(columns=mixed_column_names)
        else:
            print("no inconsistent column type found.")
            return obs

    @staticmethod
    def get_barcode2group(
        obs: pd.DataFrame,
        annot: pd.DataFrame,
        fromcol: str = "cl",
        tocol: List[str] = ["subclass_id"],
        index_col: str = "barcode",
    ):
        """Generate by GPT-4"""
        # Create a DataFrame with the mapping
        f2t = annot.loc[:, [fromcol] + tocol].copy()

        from_obs = obs[fromcol]
        if from_obs.dtype.name == "category":
            print(f"{fromcol} is category, change to int")
            obs[fromcol] = from_obs.astype(int)

        # Merge the mapping DataFrame with the original obs DataFrame
        result = obs.merge(f2t, how="left", on=fromcol)
        result.set_index(index_col, drop=False, inplace=True)
        return result

    @staticmethod
    def get_rawcount_ann(ann: ad.AnnData, outfnm: str | None = None) -> ad.AnnData:
        new_ann = ann.copy()
        # remove layers in new_ann
        new_ann.layers = None
        new_ann.X = ann.layers["rawcount"]
        new_obs = Allen.check_pandas_column_type_consistent(obs=new_ann.obs)
        new_ann.obs = new_obs
        if outfnm is not None:
            new_ann.write(outfnm)
        return new_ann

    @staticmethod
    def generate_reduced_ann(ann: ad.AnnData, outfnm: str | None = None) -> ad.AnnData:
        new_ann = ann.copy()
        new_ann.layers = None
        new_obs = Allen.check_pandas_column_type_consistent(obs=new_ann.obs)
        new_ann.obs = new_obs
        if outfnm is not None:
            new_ann.write(outfnm)
        return new_ann

def rename_allensubclass(sc: str) -> str:
    return sc.replace(" ", '_').replace("/", "-")

final_peak_bed_fnm = os.path.join(
    proj_dir, "supple.07.peakcalling.allinone",
    "mba.whole.sa2.final.peak.srt.bed")

def read_bed4(fnm) -> pd.DataFrame:
    bed: pd.DataFrame = pd.read_csv(fnm, sep = "\t",
                      header = None,
                      index_col = False)
    bed.columns = ['chrom', 'start', 'end', 'name']
    bed.set_index('name', drop = False, inplace = True)
    return bed
    
