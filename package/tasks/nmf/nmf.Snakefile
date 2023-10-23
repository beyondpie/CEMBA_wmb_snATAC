configfile : "config.yaml"
import os

n_rerun = config['n_rerun']
index_rerun = list(range(n_rerun))
max_index = max(index_rerun)
min_index = min(index_rerun)

# for test
# n_modules = [20, 40]
# modules for subclass
# n_modules = [40, 60, 80, 90, 100, 110, 120, 150]
# moduels for L3
# n_modules = [40, 60, 80, 100, 120, 150, 160, 180]

mod_from = config['mod_from']
mod_to = config['mod_to']
mod_by = config['mod_by']
n_modules = list(range(mod_from, mod_to, mod_by))

system = config["system"]
python_bin = config['python'][system]
Rscript_bin = config['Rscript'][system]
homer_bin = config['homer'][system]

local_dir = config["local_dir"]
code_dir = f'{config["code_dir"][system]}/{local_dir}'
work_dir = f'{config["work_dir"][system]}/{local_dir}'

mat_h5_file = f'{work_dir}/{config["mat_pbyc_h5"]}'
peak_nm_file = f'{work_dir}/{config["peak_nm_file"]}'
cluster_nm_file = f'{work_dir}/{config["cluster_nm_file"]}'

tag = config['tag']
out = config['out']
log_dir = f"{work_dir}/{out}/log"
flag_dir = f"{work_dir}/{out}/flag"
nmf_dir = f"{work_dir}/{out}/out"
os.makedirs(nmf_dir, exist_ok = True)

rule all:
    input:
        # expand("{d}/nmf.data.prepare.done", d = flag_dir)
        expand("{d}/nmf.{r}.{n}.nmf.done",
               d = flag_dir, r = n_modules, n = index_rerun),
        expand("{d}/nmf.{r}.post_nmf.done", d = flag_dir, r = n_modules),
        expand("{d}/nmf.{r}.stat.done", d = flag_dir, r = n_modules),
        # expand("{d}/nmf.sum.done", d = flag_dir),
        expand("{d}/nmf.plot.{r}.done", d = flag_dir, r = n_modules)

# rule prepare_data:
#     output:
#         done = touch(expand("{d}.nmf.data.prepare.done", d = flag_dir))
#     log:
#         expand("{d}/nmf.data.prepare.done", d = flag_dir)
#     shell:
#         """
#         {Rscript_bin} {code_dir}/01.prepare.cpm.R \
#             --cpm {cbyp_cpm_file} \
#             --needFilter {need_filter} \
#             --filterPeakFile {filtered_peak_file} \
#             --outPrefix {out_prefix}
#             --outDir {data_dir} 2> {log}
#         """
        
rule nmf:
    input:
        mat_h5_file
    output:
        done = touch(expand("{d}/nmf.{{r}}.{{n}}.nmf.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.{{n}}.log", d = log_dir)
    shell:
        """
        {python_bin} {code_dir}/02.nmf.py -p {input} \
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
        {python_bin} {code_dir}/02.post_nmf.py \
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
        {Rscript_bin} {code_dir}/02.nmfATAC.plotH.R \
            -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.H.mx \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        {Rscript_bin} {code_dir}/02.nmfATAC.plotW.R \
            -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.W.mx  \
           --output {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        {python_bin} {code_dir}/02.nmfATAC.stat.py -p {input.mat} \
           -x {input.peak} -y {input.cluster} \
           --basis {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.W.mx \
           --coef {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.H.mx \
            -c 0.2 \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} 2> {log}
        
        ## How can we remove the >>
        {Rscript_bin} {code_dir}/02.nmfATAC.statBox.R \
           -i {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index}.statH \
           -o {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.n{min_index} \
           >> {nmf_dir}/nmfPmat.{tag}.r{wildcards.r}.sta.txt
        """

# rule sumNMF:
#     input:
#         expand("{d}/nmf.{r}.stat.done", d = flag_dir, r = n_modules)
#     output:
#         touch(f"{flag_dir}/nmf.sum.done")
#     log:
#         f"{log_dir}/nmf.sum.log"
#     shell:
#         """
#         {Rscript_bin} {code_dir}/03.sumnmf.R \
#              --nmfDir {nmf_dir} --tag {tag}  2> {log}
#         """
        
rule plotNMF:
    input:
        expand("{d}/nmf.{{r}}.stat.done", d = flag_dir)
    output:
        touch(expand("{d}/nmf.plot.{{r}}.done", d = flag_dir))
    log:
        expand("{d}/nmf.{{r}}.plot.log", d = log_dir)
    shell:
        """
        {Rscript_bin} {code_dir}/04.nmf.plot.R --nmfDir {nmf_dir} \
            --module {wildcards.r} \
            --tag {tag} 2> {log}
        """
