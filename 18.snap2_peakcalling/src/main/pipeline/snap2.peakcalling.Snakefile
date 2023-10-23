# We have subset into files since later peakcalling needs on oasis,
# where it has slow IO. So we prepare the files before that on condo.
# then we can copy data there or just on condo.

envvars:
    "PATH"

import os
import pandas as pd
from typing import List, Dict
import pyprojroot

system = config["sys"]
proj_dir = str(pyprojroot.here())
work_dir = os.path.join(proj_dir, "18.snap2_peakcalling")
python_script_dir = os.path.join(
    proj_dir, "18.snap2_peakcalling", "src/main/python")
rdir = os.path.join(work_dir, "src/main/resource")
fdir = os.path.join(work_dir,"out", system, "flag")
ldir = os.path.join(work_dir,"out", system, "log")

## since macs2 will generate lots of files
## use oasis to store the data for tscc
out_macs2_dir = config["out_macs2"]
odir_macs2 = os.path.join(out_macs2_dir, "out", system, "macs2")

for i in [rdir, odir_macs2, fdir, ldir]:
    os.makedirs(i, exist_ok = True)

with open(f"{rdir}/mba.whole.sample.lst", 'r') as f:
    samples = [l.strip() for l in f.readlines()]
snap2s_dir = os.path.join(proj_dir, "17.snapatac2", "sa2_qc_dlt", "rm_dlt")
snap2_files = [f"{snap2s_dir}/{s}_rm_dlt.h5ad" for s in samples]

# one group could have multiple clusters
# so we call this for one group with multiple clusters.

group = config["gp"]
cluster2size_fnm = os.path.join(rdir, config["cluster2size_fnm"])
cluster2size: pd.DataFrame = pd.read_csv(
    cluster2size_fnm, sep = ',', header = None,
    names = ["cluster", "size", "early_size", "late_size"],
    index_col = None)
clusters: List[str] = cluster2size['cluster'].to_list()
cluster2size.set_index('cluster', inplace = True)

# bedtag: List[str] = ["all", "biorep", "pseudorep"]
bedfile_bedtag :str = config["bedtag"]
if not ',' in bedfile_bedtag:
    bedtag: List[str] = [bedfile_bedtag]
else:
    bedtag: List[str] = bedfile_bedtag.split(",")

all_bedtags = ["all", "biorep", "pseudorep"]
bed_outdirs: Dict[str, str] = {t: os.path.join(work_dir, "out", system, f"{t}_bed")
              for t in all_bedtags}

for i in bed_outdirs.values():
    os.makedirs(i, exist_ok = True)
def get_bed_outdir(wildcards) -> str:
    return bed_outdirs[wildcards.bedtag]

# def get_size(wildcards) -> List[int]:
#     return cluster2size.loc[wildcards.cluster,
#                             ['size', 'early_size', 'late_size']]

def get_bed_clusters(wildcards) -> List[str]:
    if wildcards.bedtag != "all":
        return ['NA']
    else:
        return clusters
    
macs2 = config["macs2"]
nmin_macs2 = config["nmin_macs2"]

debug = config["debug"]
queue = config["queue"]
walltime = int(config["walltime"])
if queue == "glean":
    if walltime > 8:
        walltime: int = 8
thread = config["thread"]

rule all:
    input:
        # expand("{fdir}/{group}_{bedtag}_bedfiles.done",
        #        fdir = fdir, group = group,
        #        bedtag = bedtag),
        expand("{fdir}/{group}_{cluster}.peakcalling.done",
               fdir = fdir, group = group, cluster = clusters)


# run all, biorep, pseudorep independently
# rule getBedFile:
#     input:
#         snap2 = f"{rdir}/cemba.snap2.with.fragment.hdf5",
#         b2c = f"{rdir}/{group}_barcode2cluster_bedtag-{{bedtag}}.csv"
#     output:
#         tag = touch(f"{fdir}/{group}_{{bedtag}}_bedfiles.done")
#     log:
#         f"{ldir}/prepare_{{bedtag}}_bedfiles.log"
#     threads: thread
#     resources:
#         walltime = walltime,
#         queue = queue
#     params:
#         seed = 0,
#         outdir = get_bed_outdir,
#         prefix = f"{group}.",
#         suffix = ".bed.gz",
#         clusters = get_bed_clusters,
#         debug = debug
#     script:
#         f"{python_script_dir}/prepare_bedfiles.py"

# put pseudo-rep and rep into this.
rule runMACS2:
    input:
        bedfile_tags = expand("{f}/{g}_{tag}_bedfiles.done",
                              f = fdir, g = group, tag = all_bedtags)
    output:
        tag = touch(f"{fdir}/{group}_{{cluster}}.peakcalling.done")
    log:
        f"{ldir}/{group}_{{cluster}}.macs2.log"
    threads: thread
    resources:
        walltime = walltime,
        queue = queue,
        mail="ae",
        email="debug.pie@gmail.com"
    params:
        debug = debug,
        bedfile_dirs = bed_outdirs,
        prefix = f"{group}.",
        outdir = odir_macs2,
        nmin = nmin_macs2,
        macs2 = macs2,
    script:
        f"{python_script_dir}/run_macs2.py"
