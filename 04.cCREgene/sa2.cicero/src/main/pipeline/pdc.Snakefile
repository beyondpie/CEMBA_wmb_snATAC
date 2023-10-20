import os
import sys
import pyprojroot

proj_dir = str(pyprojroot.here())
configfile : "config.yaml"

# Rscript_cicero = "/projects/ps-renlab2/szu/miniconda3/envs/cicero/bin/Rscript"
Rscript_cicero = "/home/szu/miniforge3/envs/r/bin/Rscript"
code_project_dir = "/projects/ps-renlab2/szu/projects/CEMBA2"
# work_project_dir = "/oasis/tscc/scratch/szu/projects/CEMBA2"
cCREgene_dir = "04.cCREgene/sa2.cicero"

meta_col = "sa2subclass"

code_dir = f'{code_project_dir}/{cCREgene_dir}/src/main'
# work_dir = f'{work_project_dir}/{cCREgene_dir}'
cicero_dir = os.path.join(
    proj_dir, "04.cCREgene", "sa2.cicero", "out/smamde_result"
)
outdir = os.path.join(
    proj_dir, "04.cCREgene", "sa2.cicero", "out/sa2pdc")
sum_outdir = os.path.join(
    proj_dir, "04.cCREgene", "sa2.cicero", "out/sa2pdcsum")

group_file = f"{code_project_dir}/meta/sa2.subclass.srt.txt"
with open(group_file) as f:
    lines = f.readlines()
    groups = [l.strip() for l in lines if len(l.strip()) > 1]
    groups = [g for g in groups if (not "Hypendymal_NN" in g)]
#groups = ["TU-ARH_Otp_Six6_Gaba"]

print(f"{len(groups)} groups will be used.")

k = config['k_cicero']
reduct = config['reduct_cicero']
pre_method = config['preprocess_cicero']
debug = config['debug']

peakannot_file = f'{code_project_dir}/{config["peakannot_file"]}'

cor_method = ["pearson"]
atac_subclass_cpm_file = f'{code_project_dir}/{config["atac_subclass_cpm_file"]}'
rdm_pdc_file = f'{code_project_dir}/{config["rdm_pdc_file"]}'
## NOTE: it's logCPM from Allen
allen_l2_cpm_file = f'{code_project_dir}/{config["allen_l2_cpm_file"]}'
pdc_file = f'{code_project_dir}/{config["pdc_file"]}'
chunk_size = config["chunk_size"]


rule all:
    input:
        # expand("flag/cicero_{g}.done", g = groups),
        # expand("flag/shufcicero_{g}.done", g = groups),
        expand("flag/filterCiceroByShuf_{g}.done", g = groups),
        expand("flag/filterDistalProximalConns_{g}.done", g = groups),
        "flag/mergeDistalProximalConns.done",
        expand("flag/cor_{m}.done", m = cor_method),
        expand("flag/rdmshuf_cor_{m}.done", m = cor_method)

# add pvalue, qvalue to conns but do not delete any
rule filterCiceroByShuf:
    output:
        touch("flag/filterCiceroByShuf_{group}.done")
    log:
        "log/filterCiceroByShuf_{group}.log"
    threads: 1
    resources:
        walltime="01:00:00",
        queue="condo"
    shell:
        """
        {Rscript_cicero} {code_dir}/R/03.filterCiceroByShuf.R \
            --metaCol {meta_col} \
            --group {wildcards.group} \
            --ciceroDir {cicero_dir} \
            --outdir {outdir} \
            --debug {debug} 2> {log}
        """
# filter conns based on FDR <= 0.001
# then get pdc without further filstering
rule filterDistalProximalConns:
    input:
        "flag/filterCiceroByShuf_{group}.done"
    output:
        touch("flag/filterDistalProximalConns_{group}.done")
    log:
        "log/filterDistalProximalConns_{group}.log"
    threads: 1
    resources:
        walltime="01:00:00",
        queue="condo"
    shell:
        """
        bash {code_dir}/shell/04.addTSSAnnot2Conns.sh \
          -c {outdir} -a {peakannot_file} \
          -g {wildcards.group} -m {meta_col} 2> {log}
        """

rule mergeDistalProximalConns:
    input:
        expand("flag/filterCiceroByShuf_{g}.done", g = groups)
    output:
        touch("flag/mergeDistalProximalConns.done")
    log:
        "log/mergeDistalProximalConns.log"
    shell:
        """
        bash {code_dir}/shell/05.mergeDistalProximalConns.sh \
          -c {outdir} -m {meta_col} 2> {log} 
        """
        
# NOTE: getCor depends on matFiles, which are need to be updated
# by the results from mergeDistalProximalConns
# We should include that as rule right after mergeDistalProximalConns
# and before getCor.
# Currently, I use supple.07.[01 and 02] to handle this.

rule getCor:
    # input:
    #     "flag/mergeDistalProximalConns.done"
    output:
        touch("flag/cor_{method}.done")
    log:
        "log/cor_{method}.log"
    shell:
        """
        {Rscript_cicero} {code_dir}/R/07.cor.scRNAseq.R \
            --matRNAFile {allen_l2_cpm_file} \
            --matATACFile {atac_subclass_cpm_file} \
            --pdcFile {pdc_file} \
            --method {wildcards.method} \
            --group {meta_col} \
            --outdir {sum_outdir} \
            --shuf 0 \
            --tag real --chunkSize {chunk_size} 2> {log}
        """

rule getRdmShufCor:
    # input:
    #     "flag/mergeDistalProximalConns.done"
    output:
        touch("flag/rdmshuf_cor_{method}.done")
    log:
        "log/rdmshuf_cor_{method}.log"
    shell:
        """
        {Rscript_cicero} {code_dir}/R/07.cor.scRNAseq.R \
            --matRNAFile {allen_l2_cpm_file} \
            --matATACFile {atac_subclass_cpm_file} \
            --pdcFile {rdm_pdc_file} \
            --method {wildcards.method} \
            --group {meta_col} \
            --outdir {sum_outdir} --shuf 1 --chunkSize {chunk_size} \
            --tag rdm 2> {log}
        """
    

