envvars:
    "PATH"

import os

with open("subclass-name.txt", 'r')  as f:
    subclasses = [l.strip() for l in f.readlines()]

rscript_bin = "/home/szu/mambaforge/bin/Rscript"
work_dir = "/projects/ps-renlab2/szu/projects/CEMBA2/04.cCREgene/sa2.cicero"
flag_dir =  os.path.join(work_dir, "flag_dir")
log_dir = os.path.join(work_dir, "log_dir")

walltime = config["walltime"]
queue = config["queue"]
if queue == "glean":
    if walltime > 8:
        walltime = 8

rule all:
    input:
        # expand(f"{f}/runcicero_{sc}.done", f = flag_dir, sc = subclasses)
        expand(f"{f}/runcicero_{sc}_shuffle.done", f = flag_dir, sc = subclasses)

rule run_cicero:
    output:
        touch(f"{flag_dir}/runcicero_{{sc}}.done")
    log:
        f"{log_dir}/runcicero_{{sc}}.log"
    threads: config["threads"]
    resources:
        walltime = walltime,
        queue = queue
    shell:
        """
        {rscript_bin} run_cicero.R {wildcards.sc} 2>&1 > {log}
        """
        
    
rule run_cicero_shuf:
    output:
        touch(f"{flag_dir}/runcicero_{{sc}}_shuffle.done")
    log:
        f"{log_dir}/runcicero_{{sc}}_shuffle.log"
    threads: config["threads"]
    resources:
        walltime = walltime,
        queue = queue
    shell:
        """
        {rscript_bin} run_cicero_shuf.R {wildcards.sc} 2>&1 > {log}
        """
