envvars:
    "PATH"

import os
import pyprojroot
import sys

proj_dir = str(pyprojroot.here())

mod = config['module']
myrange = list(range(1, mod + 1))

Rscript_bin = "/home/szu/miniforge3/envs/r/bin/Rscript"
homer_bin = "/home/szu/miniforge3/envs/r/bin/findMotifsGenome.pl"
code_project_dir = "/projects/ps-renlab2/szu/projects/CEMBA2"
motif_dir = "24.motifanalysis"

code_dir = os.path.join(code_project_dir, motif_dir, "src/main")
work_dir = os.path.join(proj_dir, "24.motifanalysis", "out")
tag = config['tag']
nmf_dir = os.path.join(proj_dir, "20.nmf", "out",
                       f"{tag}_nmf", "out")
out = config['out']
out_dir = os.path.join(work_dir, out)
log_dir = f"{work_dir}/{out}/log"
flag_dir = f"{work_dir}/{out}/flag"
motif_dir = f"{out_dir}/nmf.{tag}.r{mod}.motif"
for d in [out_dir, log_dir, flag_dir, motif_dir]:
    os.makedirs(d, exist_ok = True)
njob = config["njob"]
nmf_bed_dir = os.path.join(nmf_dir, f"nmf.{tag}.r{mod}.motif")

rule all:
    input:
        f"{flag_dir}/nmf.{tag}.splitPeakByModule.{mod}.done",
        expand("{d}/nmf.{t}.homer.{m}.n{i}.done",
               d = flag_dir, t = tag, m = mod, i = myrange)


rule splitPeakByModule:
    output:
        touch(f"{flag_dir}/nmf.{tag}.splitPeakByModule.{mod}.done")
    log:
        f"{log_dir}/nmf.{tag}.splitPeakByModule.{mod}.log"
    shell:
        """
        {Rscript_bin} {code_dir}/R/05.splitPeakByModule.R --nmfDir {nmf_dir} \
            --module {mod} --tag {tag} 2> {log}
        """
rule findMotif:
    input:
        f"{flag_dir}/nmf.{tag}.splitPeakByModule.{mod}.done"
    output:
        touch(expand("{d}/nmf.{t}.homer.{m}.n{{i}}.done", d = flag_dir, t = tag, m = mod))
    log:
        expand("{d}/nmf.{t}.homer.{m}.n{{i}}.log", d = log_dir, t = tag, m = mod)
    shell:
        """
        {homer_bin} {nmf_bed_dir}/r{mod}_n{wildcards.i}.cCREs.bed \
             mm10 {motif_dir}/homer_n{wildcards.i} \
             -nomotif -size given -p {njob} 2> {log}
        """
