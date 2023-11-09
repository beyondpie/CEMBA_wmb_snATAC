envvars:
    "PATH"
import os
import pyprojroot
import pandas as pd
from typing import List

rscript = "/home/szu/mambaforge/envs/seurat/bin/Rscript"
proj_dir = str(pyprojroot.here())
work_dir = os.path.join(proj_dir, "18.snap2_peakcalling")
script_dir = os.path.join(work_dir, "src/main")
rdir = os.path.join(work_dir, "src/main/resource")
blacklist = os.path.join(proj_dir, "meta", "mm10.blacklist.bed.gz")
chrom_size = os.path.join(proj_dir, "meta", "mm10.chrom.sizes")

oasis_dir = os.path.join("/oasis/tscc/scratch/szu/projects/CEMBA2",
                      "18.snap2_peakcalling")
macs2_dir = os.path.join(oasis_dir, "out/tscc", "macs2")
for i in [blacklist, chrom_size, macs2_dir]:
    if not os.path.exists(i):
        raise RuntimeError(f"{i} does not exist.")
    
outd_rpdpeak = os.path.join(oasis_dir, "out/tscc/rpdpeak")
outd_iter_mergepeak = os.path.join(oasis_dir, "out/tscc/itermerge")
fdir = os.path.join(oasis_dir, "out", "tscc/flag")
ldir = os.path.join(oasis_dir, "out", "tscc/log")

outd_sum_rpdpeak = os.path.join(work_dir, "out/tscc", )
outd_sum_mergepeak = os.path.join(work_dir, "out/tscc")
outd_unionpeak = os.path.join(
    work_dir,"out/tscc/unionpeak")

for d in [outd_rpdpeak, outd_iter_mergepeak,
          fdir, ldir, outd_unionpeak]:
    os.makedirs(d, exist_ok = True)

def get_pL4(L4_fnm:str, prefix:str = "nn.") -> List[str]:
    d: pd.DataFrame = pd.read_csv(
        L4_fnm, sep = ",", header = None,
        names = ['cluster', 'size', 'early_size', 'later_size'],
        index_col = None
    )
    r = [f"{prefix}{i}" for i in d['cluster'].to_list()]
    return r

nn_L4_fnm = os.path.join(rdir, "nn_L4pc2sizes_cca-k49-cl_v1.csv")
nn_pL4s: List[str] = get_pL4(nn_L4_fnm, prefix = "nn.")
neuron_L4_fnm = os.path.join(rdir, "neuron_L4pc2sizes_cca-k50-cl_v1.csv")
neuron_pL4s: List[str] = get_pL4(neuron_L4_fnm, prefix = "neuron.")
pL4s: List[str] = nn_pL4s + neuron_pL4s

print(f"{len(pL4s)} pL4 cls for merging peaks.")

rule all:
    input:
        expand("{f}/{cl}_rpdpeak.done", f = fdir, cl = pL4s),
        f"{fdir}/sum_rpdpeak.done",
        f"{fdir}/iter_merge_peak.done",
        expand("{f}/{cl}_unionpeak.done",
               f = fdir, cl = pL4s)

rule get_rpdpeak:
    output:
        touch(f"{fdir}/{{cl}}_rpdpeak.done")
    log:
        f"{ldir}/{{cl}}_rpdpeak.log"
    threads: 1
    resources:
        walltime = 1,
        queue = "glean",
        mail = "a",
        email="debug.pie@gmail.com"
    conda: "sa22"
    shell:
        """
        bash {script_dir}/shell/get_reproduce_peak_within_cluster.sh \
          -i {wildcards.cl} \
          -p {macs2_dir} \
          -o {outd_rpdpeak} 2> {log}
        """
    
rule sum_rpdpeak:
    input:
        expand("{f}/{cl}_rpdpeak.done", f = fdir, cl = pL4s)
    output:
        touch(f"{fdir}/sum_rpdpeak.done")
    log:
        f"{ldir}/sum_rpdpeak.log"
    threads: 1
    resources:
        walltime = 1,
        queue = "glean",
        mail = "a",
        email="debug.pie@gmail.com"
    conda: "seurat"
    shell:
        """
        {rscript} {script_dir}/R/sumReproducePeaks.R {outd_rpdpeak} {outd_sum_rpdpeak} \
            2> {log}
        """
rule iter_mergepeak:
    input:
        f"{fdir}/sum_rpdpeak.done"
    output:
        touch(f"{fdir}/iter_merge_peak.done")
    log:
        f"{ldir}/iter_mergepeak.log"
    threads: 5
    resources:
        walltime = 2,
        queue = "glean",
        mail = "a",
        email = "debug.pie@gmail.com"
    conda: "seurat"
    shell:
        """
        {rscript} {script_dir}/R/iterativeMergePeak.R \
           -i {outd_sum_rpdpeak}/mba.whole.naiveSummitList.list \
           -g mm10 --extend 500 \
           --blacklist  {blacklist} --chromSize {chrom_size}  \
           --ncpu {threads} -d {outd_iter_mergepeak} \
           --dsum {outd_sum_mergepeak} -o mba.whole 2> {log}
        """

rule export_unionpeak:
    input:
        f"{fdir}/iter_merge_peak.done"
    output:
        f"{outd_sum_mergepeak}/mba.whole.union.peak.bed"
    threads: 1
    resources:
        walltime = 1,
        queue = "glean"
    shell:
        """
        bash {script_dir}/shell/export_unionpeak.sh \
           {outd_sum_mergepeak}/mba.whole.filteredNfixed.union.peakSet \
           {output}
        """
        
rule intersect_mergepeak:
    input:
        f"{outd_sum_mergepeak}/mba.whole.union.peak.bed"
    output:
        touch(f"{fdir}/{{cl}}_unionpeak.done")
    threads: 1
    resources:
        walltime = 1,
        queue = "glean",
        mail = "a",
        email = "debug.pie@gmail.com"
    conda: "sa22"
    shell:
        """
        bash {script_dir}/shell/intersect_mergepeak.sh \
           {input} \
           {outd_iter_mergepeak}/{wildcards.cl}.filterNfixed.peakset \
           {outd_unionpeak}/{wildcards.cl}.union.peak.bed
        """
