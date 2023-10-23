#envvars:
#    "PATH"
import os
import pyprojroot
import pandas as pd
from typing import List

proj_dir = str(pyprojroot.here())
work_dir = os.path.join(proj_dir, "18.snap2_peakcalling")
script_dir = os.path.join(work_dir, "src/main")

union_bed_file = os.path.join(work_dir, "out/tscc",
                               "mba.whole.union.peak.srt.bed")
if not os.path.exists(union_bed_file):
    raise FileNotFoundError(union_bed_file)

rnd_bed_file = os.path.join(work_dir, "out/tscc",
                               "mba.whole.shuffle.removeovlp.bed")
if not os.path.exists(rnd_bed_file):
    raise FileNotFoundError(rnd_bed_file)

with open(f"{proj_dir}/17.snapatac2/meta/mba.whole.sample.lst", 'r') as f:
    samples = [l.strip() for l in f.readlines()]

snap2_dir = os.path.join(proj_dir,
                         "17.snapatac2",
                         "sa2_qc_dlt/rm_dlt")

fdir = os.path.join(work_dir, "out", "tscc/sa2pmat" , "flag")
ldir = os.path.join(work_dir, "out", "tscc/sa2pmat" , "log")
odir_union = os.path.join(work_dir, "out", "tscc/sa2pmat" , "union_pmat")
odir_rnd = os.path.join(work_dir, "out", "tscc/sa2pmat" , "union_pmat_rnd")

for d in [fdir, ldir, odir_union, odir_rnd]:
    os.makedirs(d, exist_ok = True)

rule all:
    input:
        expand("{f}/{s}.sa2pmat.flag", f=fdir, s=samples),
        expand("{f}/{s}.sa2pmat_rnd.flag", f=fdir, s=samples)

rule sa2pmat_union:
    input:
        bedfnm=union_bed_file
    output:
        touch(f"{fdir}/{{s}}.sa2pmat.flag")
    log:
        f"{ldir}/{{s}}.sa2pmat.log"
    params:
        snap2_dir = snap2_dir,
        suffix = "_rm_dlt.h5ad",
        out_dir = odir_union,
        out_suffix = "_union_pmat.h5ad"
    #conda: "sa22"
    threads: 2
    resources:
        walltime="1:00:00",
        queue="glean"
    script:
        f"{script_dir}/python/sa2pmat.py"
        
rule sa2pmat_rnd:
    input:
        bedfnm=rnd_bed_file
    output:
        touch(f"{fdir}/{{s}}.sa2pmat_rnd.flag")
    log:
        f"{ldir}/{{s}}.sa2pmat_rnd.log"
    params:
        snap2_dir = snap2_dir,
        suffix = "_rm_dlt.h5ad",
        out_dir = odir_rnd,
        out_suffix = "_rnd_pmat.h5ad"
    #conda: "sa22"
    threads: 2
    resources:
        walltime="1:00:00",
        queue="glean"
    script:
        f"{script_dir}/python/sa2pmat.py"
