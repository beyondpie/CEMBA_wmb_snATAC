configfile: "config.yaml"
import os
import sys
import pyprojroot

proj_dir = str(pyprojroot.here())

system = config['system']
python_bin = config["python_bin"][system]
code_project_dir = config['code_project_dir'][system]
work_project_dir = config['work_project_dir'][system]
local_dir = config['local_dir']
code_dir = f'{code_project_dir}/{local_dir}/src/main'

pdc_dir = os.path.join(proj_dir, "04.cCREgene", "sa2.cicero",
                       "out", "sa2pdc")
group_file = os.path.join(proj_dir, "meta", "sa2.subclass.srt.txt")
with open(group_file) as f:
    lines = f.readlines()
    groups = [l.strip() for l in lines if len(l.strip()) > 1]
    groups = [g for g in groups if (not "Hypendymal_NN" in g)]
# groups = ["TU-ARH_Otp_Six6_Gaba"]
print(f"{len(groups)} groups will be used.")


tfscan_dir = f'{work_project_dir}/{local_dir}/{config["tfscan_dir"]}'
os.makedirs(tfscan_dir, exist_ok = True)

# pdc_suffix = f'{code_project_dir}/{local_dir}/{config["pdc_suffix"]}'
# with open(pdc_suffix, 'r') as f:
#     suffix = [l.strip() for l in f.readlines() if len(l.strip()) > 1]

threshold_gimme = config["threshold_gimme"]

# subclass4GRN_file = f'{code_project_dir}/{local_dir}/{config["subclass4GRN"]}'
# with open(subclass4GRN_file, 'r') as f:
#     subclases = [l.strip() for l in f.readlines() if len(l) > 1]

scRNA_dir = f'{work_project_dir}/{local_dir}/{config["allenRNA_dir"]}'
GRN_dir = f'{work_project_dir}/{local_dir}/{config["GRN_dir"]}'
GRN_all = os.path.join(GRN_dir, "baseGRN_all")
GRN_sc = os.path.join(GRN_dir, "baseGRN_subclass")
for d in [GRN_dir, GRN_all, GRN_sc]:
    os.makedirs(d, exist_ok = True)

njobs=config['njobs']


rule all:
    input:
        expand("flag/gimmemotifs.{g}.done", g = groups),
        expand("flag/runGRN.{g}.baseGRNall.done", g = groups),
        # expand("flag/runGRN.{g}.baseGRNsc.done", g = groups)

# install mm10 genome firstly following supple script.
rule runGimmemotifs:
    output:
        touch("flag/gimmemotifs.{group}.done")
    log:
        "log/gimmemotifs.{group}.log"
    threads: 1
    shell:
        """
        {python_bin} {code_dir}/python/01.runGimmemotifs.py \
            --pdcbedpe \
             {pdc_dir}/sa2subclass.{wildcards.group}.pdc.bedpe \
            --outdir {tfscan_dir} \
            --threshold {threshold_gimme} 2> {log}
        """


# use 02.mergeGimme to merge all the results into one.
# use 03.seurat2anndata.py to generate the scRNAseq for GRN.


rule runGRN_baseGRNall:
    output:
        touch("flag/runGRN.{group}.baseGRNall.done")
    threads: njobs
    log:
        "log/runGRN.{group}.baseGRNall.log"
    shell:
        """
        {python_bin} {code_dir}/python/04.runGRN.py --tfdir {tfscan_dir} \
           --scRNAdir {scRNA_dir} --outdir {GRN_all} \
           --group {wildcards.group} --njobs {njobs} \
           --baseGRN_all 1 > {log}
        """


rule runGRN_baseGRNsc:
    output:
        touch("flag/runGRN.{group}.baseGRNsc.done")
    threads: njobs
    log:
        "log/runGRN.{group}.baseGRNsc.log"
    shell:
        """
        {python_bin} {code_dir}/python/04.runGRN.py --tfdir {tfscan_dir} \
           --scRNAdir {scRNA_dir} --outdir {GRN_sc} \
           --group {wildcards.group} --njobs {njobs} \
           --baseGRN_all 0 > {log}
        """
