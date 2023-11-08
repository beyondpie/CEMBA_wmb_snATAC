import os
from typing import Dict
system: str = "encoder"
project_dict: Dict[str, str] = {
    "imac": "/Users/szu/git-recipes/mouseBrainAtlas",
    "encoder": "/projects/ps-renlab/szu/projects/CEMBA2"
}
genome = "sa2default"

project_dir = project_dict[system]
rm_dlt_dir = f"{project_dir}/17.snapatac2/sa2_qc_dlt/rm_dlt"
with open(f"{project_dir}/17.snapatac2/meta/mba.whole.sample.lst", 'r') as f:
    samples = [l.strip() for l in f.readlines()]
# test only
# samples = ["CEMBA171206_3C", "CEMBA171207_3C"]
# samples = ["CEMBA171206_3C"]

out_dir = f"{project_dir}/17.snapatac2/sa2_{genome}_gmat"
log_dir = f"{out_dir}/log"
flag_dir = f"{out_dir}/flag"
# sample-level gmat
sgmat_dir = f"{out_dir}/sgmat"
for d in [out_dir, log_dir, flag_dir, sgmat_dir]:
    os.makedirs(d, exist_ok = True)


def get_sample(wildcards):
    return wildcards.s

rule all:
    input:
        expand("{f}/{s}_{g}_gmat.done",
               f = flag_dir, g = genome, s = samples),
        f"{flag_dir}/{genome}_gmat_merged.done"

rule sgmat:
    input:
        snap_file = expand("{i}/{{s}}_rm_dlt.h5ad", i = rm_dlt_dir)
    output:
        gmat_file = expand("{o}/{{s}}_{g}_gmat.h5ad",
                           o = sgmat_dir, g = genome),
        tag = touch(expand("{f}/{{s}}_{g}_gmat.done",
                           f = flag_dir, g = genome))
    log:
        expand("{l}/{{s}}_{g}_gmat.log", l = log_dir, g = genome)
    params:
        sample = get_sample,
        genome = genome
    threads: 1
    script:
        f"{project_dir}/17.snapatac2/script/sa2.get.sample.gmat.py"

rule merge_sgmat:
    input:
        snap_files = expand("{o}/{s}_{g}_gmat.h5ad",
                            o = sgmat_dir, s = samples, g = genome)
    output:
        merge_snap = f"{out_dir}/{genome}_gmat_merged.h5ad",
        tag = touch(f"{flag_dir}/{genome}_gmat_merged.done")
    params:
        genome = genome
    log:
        f"{log_dir}/{genome}_gmat_merged.log"
    threads: 4
    script:
        f"{project_dir}/17.snapatac2/script/sa2.merge.gmat.py"
    
        
