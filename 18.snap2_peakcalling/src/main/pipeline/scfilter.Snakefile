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
odir = os.path.join(work_dir, "out/scfilter")
fdir = os.path.join(odir, "flag")
ldir = os.path.join(odir, "log")
odir_peakfrac = os.path.join(odir, "peakfrac")
odir_fitbgmodel = os.path.join(odir, "fitfrac_bg")
odir_clfilter = os.path.join(odir, "clfileter")

for d in [odir, fdir, ldir, odir_peakfrac, odir_fitbgmodel,
          odir_clfilter]:
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
print(f"{len(pL4s)} pL4 cls for filtering peaks.")

# pL4ids = list(range(1,3))
pL4ids = list(range(1,1464))

# sa2pmatd = os.path.join(work_dir, "out/tscc/sa2pmat")
# sa2pmatd_union = os.path.join(sa2pmatd, "union_pmat")
# sa2pmatd_rnd = os.path.join(sa2pmatd, "union_pmat_rnd")

rule all:
    input:
        # f"{fdir}/peakfrac.done",
        expand("{f}/fitbgmodel_{cl}.done", f = fdir, cl = pL4ids)

# let's do this in our script interactively.
# rule get_peakfrac:
#     output:
#         touch(f"{fdir}/peakfrac.done")
#     params:
#         sa2pmatd_fg = sa2pmatd_union,
#         sa2pmatd_bg = sa2pmatd_rnd,
#         outdir = odir_peakfrac
#     threads:
#         20
#     resources:
#         walltime = 8,
#         queue = "glean",
#         mail = "a",
#         tag = "icelake:mem1024"
#     conda:
#         "sa22"
#     script:
#         "{script_dir}/python/sa2_get_peakfrac.py"

rule fit_bgmodel:
    # input:
    #     f"{fdir}/peakfrac.done"
    output:
        touch(f"{fdir}/fitbgmodel_{{cl}}.done")
    log:
        f"{ldir}/fitbgmodel_{{cl}}.done"
    threads: 5
    resources:
        walltime = 4,
        queue = "condo",
        mail = "a"
    conda:
        "seurat"
    shell:
        """
        {rscript} {script_dir}/R/fitbgmodel.R {wildcards.cl} 2> {log}
        """
        
        
        
