import os
a = ["a1", "a2", "a3"]

out_dir = "out"
log_dir = "log"
for i in [out_dir, log_dir]:
    os.makedirs(i, exist_ok = True)

script_dir = "../R"

rule all:
    input:
        expand("{o}/{t}.out", o = out_dir, t = a)

rule getOut:
    input:
        # afnm = lambda w: f"{out_dir}/{w.t}.in",
        # afnm = expand("{o}/{{t}}.in", o = out_dir),
        afnm = f"{out_dir}/{{t}}.in",
        # bfnm = lambda w: f"{out_dir}/{w.t}.in"
        bfnm = f"{out_dir}/{{t}}.in"
        # bfnm = expand("{o}/{{t}}.in", o = out_dir)
    output:
        touch(expand("{o}/{{t}}.out", o = out_dir))
    log:
        f"{log_dir}/{{t}}.log"
        # wildcard in log filed does not work
        # fnm = lambda w: f"{log_dir}/{w.t}.log"
    params:
        a = [1,2,3]
    script:
        f"{script_dir}/test.snakemake.wildcards.R"
        
        
