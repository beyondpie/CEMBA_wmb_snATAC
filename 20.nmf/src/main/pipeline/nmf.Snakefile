configfile : "config.yaml"
envvars:
    "PATH"
    
import os
import pyprojroot

n_rerun = config['n_rerun']
index_rerun = list(range(n_rerun))
max_index = max(index_rerun)
min_index = min(index_rerun)

mod_from = config['mod_from']
mod_to = config['mod_to']
mod_by = config['mod_by']
n_modules = list(range(mod_from, mod_to, mod_by))

system = config["system"]
python_bin = config['python'][system]
Rscript_bin = config['Rscript'][system]
homer_bin = config['homer'][system]

code_dir = os.path.join(config["code_dir"][system],
                        "20.nmf", "src/main")
work_dir = os.path.join(config["work_dir"][system],
                        "20.nmf")

mat_h5_file = os.path.join(work_dir, config["mat_pbyc_h5"])
peak_nm_file = os.path.join(work_dir, config["peak_nm_file"])
cluster_nm_file = os.path.join(work_dir, config["cluster_nm_file"])
subclass_order_meta = os.path.join(
    config["code_dir"][system], config["subclass_order_meta"])

for f in [mat_h5_file, peak_nm_file, cluster_nm_file]:
    if not os.path.exists(f):
        raise FileNotFoundError(f)

tag = config['tag']
out = config['out']
log_dir = f"{work_dir}/{out}/log"
flag_dir = f"{work_dir}/{out}/flag"
nmf_dir = f"{work_dir}/{out}/out"

for d in [log_dir, flag_dir, nmf_dir]:
    os.makedirs(d, exist_ok = True)


rule all:
    input:
        expand("{d}/nmf.{r}.{n}.nmf.done",
               d = flag_dir, r = n_modules, n = index_rerun),
        expand("{d}/nmf.{r}.post_nmf.done",
               d = flag_dir, r = n_modules),
        expand("{d}/nmf.{r}.stat.done", d = flag_dir, r = n_modules),
        expand("{d}/nmf.sum.done", d = flag_dir),
        expand("{d}/nmf.plot.{r}.done", d = flag_dir, r = n_modules)
        
rule nmf:
    input:
        mat_h5_file
    output:
        touch(expand("{d}/nmf.{{r}}.{{n}}.nmf.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.{{n}}.log", d = log_dir)
    shell:
        """
        {python_bin} {code_dir}/python/02.nmf.py -p {input} \
            -r {wildcards.r} -n {wildcards.n} \
            -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{wildcards.n} \
            2> {log}
        """

rule post_nmf:
    input:
        expand("{d}/nmf.{r}.{n}.nmf.done",
               d = flag_dir, r = n_modules, n = index_rerun)
    output:
        touch(expand("{d}/nmf.{{r}}.post_nmf.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.post_nmf.log", d = log_dir)
    shell:
        """
        {python_bin} {code_dir}/python/02.post_nmf.py \
           -r {wildcards.r} \
           -i {nmf_dir}/nmfPmat.{tag} \
           -o {nmf_dir}/nmfPmat.{tag} \
           -d {min_index} -s {min_index} -e {max_index} 2> {log}
        """
        
rule stat:
    input:
        depend = rules.post_nmf.output,
        peak = peak_nm_file,
        cluster = cluster_nm_file,
        mat = mat_h5_file
    output:
        done = touch(expand("{d}/nmf.{{r}}.stat.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.stat.log", d = log_dir)
    shell:
        """
        {Rscript_bin} {code_dir}/R/02.nmfATAC.plotH.R \
            -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.H.mx \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        {Rscript_bin} {code_dir}/R/02.nmfATAC.plotW.R \
            -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.W.mx  \
           --output {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        {python_bin} {code_dir}/python/02.nmfATAC.stat.py -p {input.mat} \
           -x {input.peak} -y {input.cluster} \
           --basis {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.W.mx \
           --coef {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.H.mx \
            -c 0.2 \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        ## How can we remove the >>
        {Rscript_bin} {code_dir}/R/02.nmfATAC.statBox.R \
           -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.statH \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} \
           >> {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.sta.txt
        """

rule sumNMF:
    input:
        expand("{d}/nmf.{r}.stat.done", d = flag_dir, r = n_modules)
    output:
        touch(f"{flag_dir}/nmf.sum.done")
    log:
        f"{log_dir}/nmf.sum.log"
    shell:
        """
        {Rscript_bin} {code_dir}/R/03.sumnmf.R \
             --nmfDir {nmf_dir} --tag {tag}  2> {log}
        """
        
rule plotNMF:
    input:
        expand("{d}/nmf.{{r}}.stat.done", d = flag_dir)
    output:
        touch(expand("{d}/nmf.plot.{{r}}.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.plot.log", d = log_dir)
    shell:
        """
        {Rscript_bin} {code_dir}/R/04.nmf.plot.R --nmfDir {nmf_dir} \
            --module {wildcards.r} \
            --tag {tag} 2> {log}
        """
