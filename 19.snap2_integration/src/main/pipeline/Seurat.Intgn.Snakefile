# let snaekmake pass this env
# so that it can find the Rscript
# and qsub needs -V param (our pbs profile supports it)
envvars:
    "PATH"
    
import os
import pyprojroot

# currently, we don't need this
# config: "config.yaml"
py = "/home/szu/mambaforge/envs/sa2/bin/python"

proj_dir = str(pyprojroot.here())
rsc_dir = os.path.join(proj_dir,
                       "19.snap2_integration",
                       "src/main/resource")

# @preprocess
# NOTE: gp like: neuron, nn, or somespecific one
# we prepare the the corresponding outside of the pipeline
# - ann data
# - obs_meta as csv file
gp: str = config["gp"]

odir = os.path.join(proj_dir, "19.snap2_integration",
                       f"out/{gp}")

tf_data_dir = os.path.join(proj_dir, "19.snap2_integration",
                      f"out/transferLabel")
tfodir = os.path.join(proj_dir, "19.snap2_integration",
                      f"out/tf/{gp}")


# intgn_methods = ["rpca", "cca", "harmony", "mnn"]
# intgn_methods = ["rpca", "cca", "mnn"]
intgn_methods = config["intgn_method"].split(",")
k_anchor = [5, 30]

# @preprocess
# if group-specific features are used, we may need to
# add them to the feature list
# k8, TF, markers, vf
# vf: use variable features in the data
features = config['feature'].split(",")

n_pca = 50
# index start from 1 as in R
# use this to check if pca
# from 2:n_pca would be better
start_pca_atac = 1
start_pca_rna = 1
# allen means: 10xv3 + 10xv2 + multiome
# allen_techs = ["10xv3", "allen"]
allen_techs = config["allen_techs"].split(",")
useAllenAsRef = config["allen_as_ref"]
dsSeurat = config["ds_seu"]
## downsample lower and upper bound 
dsMin = config["ds_nmin"]
dsMax = config["ds_nmax"]
# resource
queue = "hotel"

# DEBUG config
if config["debug"] > 0:
    print("Under DEBUG Mode: ... ")
    features = ["TF"]
    allen_techs = ["10xv3"]
    intgn_methods = ["rpca"]
    k_anchors = [5]
    odir = os.path.join(proj_dir, "19.snap2_integration",
                        f"out/test/{gp}")
    rsc_dir = os.path.join(proj_dir,
                        "19.snap2_integration",
                        "src/test/resource")
else:
    print("Under Normal Mode: ...")

flag_dir = os.path.join(odir, "flag")
ldir = os.path.join(odir, "log")
# make sure out dirs exist.
for i in [odir, flag_dir, ldir]:
    os.makedirs(i, exist_ok = True)
    
allen_ann_dir = os.path.join(rsc_dir, "norawdata_allen")
atac_ann_dir = os.path.join(rsc_dir, "norawdata_atac")

all_techs = ["atac"] + allen_techs
modality = ["atac", "allen"]
# use for subset sample, we remove this as preprocess
# atac_cluster_level = config["atac_cluster_level"]
# allen_cluster_level = config["allen_cluster_level"]
script_dir = os.path.join(proj_dir, "19.snap2_integration",
                          "src/main/R")

def get_intgn_ppn(wildcards) -> int:
    if wildcards.m == "cca":
        return 32
    else:
        return 8

def get_allen_tech(wildcards) -> str:
    return wildcards.tech

def get_annfnm(wildcards) -> str:
    if wildcards.tech == "atac":
        return f"{atac_ann_dir}/{gp}_gmat_atac_ann.h5ad"
    else:
        return f"{allen_ann_dir}/{gp}_male_{wildcards.tech}_ann_noraw.h5ad"

def get_modality(wildcards) -> str:
    if wildcards.tech == "atac":
        return "atac"
    else:
        return "rna"

def get_fromPCA(wildcards) -> int:
    if wildcards.tech == "atac":
        return start_pca_atac
    else:
        return start_pca_rna

rule all:
    input:
        expand("{o}/{g}_{tech}.s5.rds",
               o = odir, g = gp, tech = all_techs),
        expand("{o}/{g}_{tech}.{ft}.s5.rds",
               o = odir, g = gp, tech = all_techs, ft = features),
        expand("{o}/{g}_{at}_{m}_{ft}_{k}.s5.rds",
               o = odir, g = gp, at = allen_techs,
               m = intgn_methods, ft = features, k = k_anchor),
        expand("{o}/{g}_{at}_{m}_{ft}_{k}.s5.umap.rds",
               o = odir, g = gp, at = allen_techs,
               m = intgn_methods, ft = features, k = k_anchor)

rule transformAnnDatatoS5:
    input:
        get_annfnm
    output:
        f"{odir}/{gp}_{{tech}}.s5.rds",
        touch(f"{flag_dir}/{gp}_{{tech}}.s5.done")
    log:
        f"{ldir}/{gp}_{{tech}}.s5.log"
    threads: 2
    params:
        modality = get_modality
    resources:
        walltime = "02:00:00",
        queue = queue
    script:
        # TODO: need to use the lastest fns in package
        f"{script_dir}/annToS5.R"

# TODO: remove this, and merge into IntegrationData
rule runPCA:
    input:
        f"{odir}/{gp}_{{tech}}.s5.rds"
    output:
        f"{odir}/{gp}_{{tech}}.{{ft}}.s5.rds",
        touch(f"{flag_dir}/{gp}_{{tech}}.{{ft}}.s5.done")
    log:
        f"{ldir}/{gp}_{{tech}}.{{ft}}.s5.log"
    threads: 6
    params:
        nPCA = n_pca,
        downSample = 0
    resources:
        walltime = "08:00:00",
        queue = queue
    script:
        f"{script_dir}/runPCA.R"

# will read wildcards m, k, ft
# IntegrationData for UMAP and consensus matrix
rule IntegrationData:
    input:
        atacS5 = f"{odir}/{gp}_atac.{{ft}}.s5.rds",
        allenS5 = f"{odir}/{gp}_{{at}}.{{ft}}.s5.rds"
    output:
        f"{odir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.rds",
        touch(f"{flag_dir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.done")
    log:
        f"{ldir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.log"
    params:
        nPCA = n_pca,
        fromATAC = start_pca_atac,
        fromRNA = start_pca_rna
    threads: get_intgn_ppn
    resources:
        walltime = "20:00:00",
        queue = queue
    script:
        f"{script_dir}/S5Int.R"

# TODO: run Leiden + silhouette
# then calculate consensus matrix in R
rule umap:
    input:
        f"{odir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.rds"
    output:
        f"{odir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.umap.rds",
        touch(f"{flag_dir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.umap.done")
    log:
        f"{ldir}/{gp}_{{at}}_{{m}}_{{ft}}_{{k}}.s5.umap.log"
    threads: 8
    params:
        nPCA = n_pca
    resources:
        walltime = "10:00:00",
        queue = queue
    script:
        f"{script_dir}/getIntUMAP.R"

