# let snaekmake pass this env
# so that it can find the Rscript
# and qsub needs -V param (our pbs profile supports it)
envvars:
    "PATH"
    
import os
import pyprojroot

proj_dir = str(pyprojroot.here())
## encoder or tscc
system = config["sys"]
rsc_dir = os.path.join(proj_dir,
                       "19.snap2_integration",
                       f"out/transferLabel_{system}")

# @preprocess
# NOTE: gp like: neuron, nn, or somespecific one
# we prepare the the corresponding outside of the pipeline
# - ann data
# - obs_meta as csv file
gp: str = config["gp"]

odir = os.path.join(proj_dir, "19.snap2_integration",
                       f"out/tf_{system}/{gp}")

# intgn_methods = ["rpca", "cca", "harmony", "mnn"]
# intgn_methods = ["rpca", "cca", "mnn"]
if not "," in config["intgn_method"]:
    intgn_methods = [config["intgn_method"]]
else:
    intgn_methods = config["intgn_method"].split(",")

if isinstance(config["k_anchors"], int):
    print("k_anchors is an integer.")
    k_anchor = [config["k_anchors"]]
elif not "," in config["k_anchors"]:
    k_anchor = [int(config["k_anchors"])]
else:
    k_anchor = [int(i) for i in config["k_anchors"].split(",")]

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
if not "," in config["allen_techs"]:
    allen_techs = [config["allen_techs"]]
else:
    allen_techs = config["allen_techs"].split(",")
    
useAllenAsRef = config["allen_as_ref"]
dsSeurat = config["ds_seu"]
## downsample lower and upper bound 
allen_dsMin = config["allen_ds_min"]
allen_dsMax = config["allen_ds_max"]
atac_dsMin = config["atac_ds_min"]
atac_dsMax = config["atac_ds_max"]
# resource
queue = config["queue"]

# DEBUG config
if config["debug"] > 0:
    print("Under DEBUG Mode: ... ")
    allen_techs = ["10xv3"]
    k_anchor = [5]
    odir = os.path.join(proj_dir, "19.snap2_integration",
                        f"out/test/tf/{gp}")
else:
    print("Under Normal Mode: ...")

flag_dir = os.path.join(odir, "flag")
ldir = os.path.join(odir, "log")
# make sure out dirs exist.
for i in [odir, flag_dir, ldir]:
    os.makedirs(i, exist_ok = True)
    
# use for subset sample, we remove this as preprocess
# atac_cluster_level = config["atac_cluster_level"]
# allen_cluster_level = config["allen_cluster_level"]
script_dir = os.path.join(proj_dir, "19.snap2_integration",
                          "src/main/R")

# CCA for 360K cells: 330G
# with about 8 hours for anchors in 1 hour

def get_intgn_ppn(wildcards) -> int:
    if wildcards.m == "cca":
        # increase cores for cca
        return 60
    else:
        return 8

walltime = "15:00:00"
if queue in ["glean", "condo"]:
    walltime = "08:00:00"

rule all:
    input:
        expand("{f}/{g}_atac.{at}_{m}_{ft}_{k}.tf.done",
               f = flag_dir, g = gp, at = allen_techs,
               m = intgn_methods, ft = features, k = k_anchor)

rule TransferLabel:
    input:
        atacSeu = f"{rsc_dir}/{gp}_atac_noraw_seurat.rds",
        allenSeu = f"{rsc_dir}/{gp}_{{at}}_noraw_seurat.rds"
    output:
        anchor = f"{odir}/{gp}_atac.{{at}}_{{m}}_{{ft}}_{{k}}.tf.anchor.rds",
        transferLabel = f"{odir}/{gp}_atac.{{at}}_{{m}}_{{ft}}_{{k}}.tf.label.rds",
        tag = touch(f"{flag_dir}/{gp}_atac.{{at}}_{{m}}_{{ft}}_{{k}}.tf.done")
    params:
        nPCA = n_pca,
        useAllenAsRef = useAllenAsRef,
        nFeatures = 2000,
        vfon = "ref",
        dsSeurat = dsSeurat,
        atac_dsMin = atac_dsMin,
        atac_dsMax = atac_dsMax,
        allen_dsMin = allen_dsMin,
        allen_dsMax = allen_dsMax,
        atacGroupBy = "L4",
        allenGroupBy = "cl"
    log:
        f"{ldir}/{gp}_atac.{{at}}_{{m}}_{{ft}}_{{k}}.tf.log"
    threads: get_intgn_ppn
    resources:
        walltime = walltime,
        queue = queue,
        tag = "icelake:mem1024"
    script:
        f"{script_dir}/TransferLabel.R"
