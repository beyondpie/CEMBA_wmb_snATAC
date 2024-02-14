import os
configfile: "snapatac2.qc.config.yaml"

system = config["system"]
conda_env = f'{config["conda"][system]}/envs/{config["conda_env"][system]}'
python = f"{conda_env}/bin/python"
project_dir = config["project_dir"][system]
work_dir = config["work_dir"][system]
out_dir = f'{work_dir}/{config["out_dir"]}'
os.makedirs(out_dir, exist_ok = True)
# prepare out_sub_dirs
os.makedirs(f"{out_dir}/raw_fragment", exist_ok = True)
os.makedirs(f"{out_dir}/raw_stat", exist_ok = True)
os.makedirs(f"{out_dir}/qc", exist_ok = True)
os.makedirs(f"{out_dir}/doublet", exist_ok = True)
os.makedirs(f"{out_dir}/bmat", exist_ok = True)
os.makedirs(f"{out_dir}/gmat", exist_ok = True)

sample2bamfile = f'{project_dir}/{config["sample2bamfile"][system]}'
with open(sample2bamfile, 'r') as f:
    lines = [l.strip() for l in f.readlines()]
    sample2bam = {l.split(',')[0]:l.split(',')[1] for l in lines}
    
samples_file = f'{project_dir}/{config["samples_file"][system]}'
with open(samples_file, 'r') as f:
    samples = [l.strip() for l in f.readlines() if len(l) > 1]
    

chrom_size_file = f'{project_dir}/{config["chrom_size_file"]}'
gtf_file = f'{project_dir}/{config["gtf_file"]}'
blacklist_file = f'{project_dir}/{config["blacklist_file"]}'

nc = config["embed_ncomps"]
ss = config["embed_sample_size"]
dm = config["distance_metric"]
nf = config["embed_nfeature"]
embed_suffix = f"{nf}_{nc}_{ss}_{dm}"
outf_l1_embed = f"{out_dir}/l1_{system}/sa2_l1embed_{embed_suffix}.hdf5"
knn = config["knn"]
km = config["knn_method"]
knn_suffix = f"{embed_suffix}_{km}_{knn}"
outf_l1_knn = f"{out_dir}/l1_{system}/sa2_l1knn_{knn_suffix}.hdf5"
debug = config["debug"]
            
rule all:
   input:
       # snapatac2_pp_done = expand("flag/{s}_sa2.pp.done",
       #                            s = samples),
       snaptac2_l1embed_done = f"flag/sa2_l1embed_{system}_{embed_suffix}.done",
       snapatac2_l1knn_done = f"flag/sa2_l1knn_{system}_km{km}_k{knn}.done"

rule snapatac2_pp:
    input:
        lambda wildcards: sample2bam[wildcards.sample]
    output:
        touch("flag/{sample}_sa2.pp.done")
    log:
        "log/{sample}_snapatac2_pp.log"
    shell:
        """
        {python} {work_dir}/sa2.01.preprocess.py --sample {wildcards.sample} \
             --outdir {out_dir} --bamfile {input} \
             --chromfile {chrom_size_file} \
             --blacklist {blacklist_file} \
             --gtfile {gtf_file} --logfile {log} --debug {debug}
        """

rule snapatac2_l1_embed:
    # input:
        # expand("flag/{s}_sa2.pp.done", s = samples)
    output:
        touch(f"flag/sa2_l1embed_{system}_{embed_suffix}.done")
    benchmark:
        f"benchmarks/sa2_l1embed_{system}_{embed_suffix}.tsv"
    log:
        f"log/sa2_l1embed_{system}_{embed_suffix}.log"
    shell:
        """
        {python} {work_dir}/sa2.03.l1_embed.py \
           --bmatdir {out_dir}/bmat \
           --nfeature {nf} \
           --outf {outf_l1_embed} \
           --samples {samples_file} \
           --blacklist {blacklist_file} \
           --samplesize {ss} \
           --ncomps {nc} \
           --dmetric {dm} \
           --logfile {log} --debug {debug} 2> {log}
        """

rule snapatac2_l1_knn:
    input:
        f"flag/sa2_l1embed_{system}_{embed_suffix}.done"
    output:
        touch(f"flag/sa2_l1knn_{system}_km{km}_k{knn}.done")
    benchmark:
        f"benchmarks/sa2_l1knn_{system}_km{km}_k{knn}.tsv"
    log:
        f"log/sa2_l1knn_{system}_km{km}_k{knn}.log"
    shell:
        """
        {python} {work_dir}/sa2.04.l1_knn.py \
           --embed_file {outf_l1_embed} \
           --knn {knn} \
           --kmethod {km} \
           --outf {outf_l1_knn} \
           --logfile {log} --debug {debug} 2> {log}
        """
