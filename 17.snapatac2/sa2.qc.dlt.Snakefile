import os
code_dir = "/projects/ps-renlab/szu/projects/CEMBA2"
with open(f"{code_dir}/17.snapatac2/meta/mba.whole.sample.lst", 'r') as f:
    samples = [l.strip() for l in f.readlines()]
# out_dir = "/oasis/tscc/scratch/szu/projects/CEMBA2/17.snapatac2/sa2_qc_dlt"
out_dir = "/projects/ps-renlab/szu/projects/CEMBA2/17.snapatac2/sa2_qc_dlt"

# test
# code_dir = "/Users/szu/git-recipes/mouseBrainAtlas/CEMBA2"
# with open(f"{code_dir}/17.snapatac2/meta/mba.test.sample", 'r') as f:
    # samples = [l.strip() for l in f.readlines()]
# out_dir = "/Users/szu/test"

log_dir = f"{out_dir}/log"
flag_dir = f"{out_dir}/flag"
os.makedirs(log_dir, exist_ok = True)
os.makedirs(flag_dir, exist_ok = True)
# frag_dir = f"{out_dir}/frag"
# frag_sum_dir = f"{out_dir}/frag_sum"
# os.makedirs(frag_dir, exist_ok = True)
# os.makedirs(frag_sum_dir, exist_ok = True)

qc_dlt_dir = f"{out_dir}/qc_dlt"
os.makedirs(qc_dlt_dir, exist_ok = True)

rm_dlt_dir = f"{out_dir}/rm_dlt"
os.makedirs(rm_dlt_dir, exist_ok = True)

def get_sample(wildcards):
    return wildcards.s

rule all:
    input:
        # expand("{f}/{s}_qc_dlt.done", f = flag_dir, s = samples)
        f"{flag_dir}/merge_cemba_all.done"

# rule qc_dlt:
#     input:
#         sample2bam_file = f"{code_dir}/meta/sample2rawbam.csv",
#         chromsizes_file = f"{code_dir}/meta/mm10.chrom.sizes.lite",
#         blacklist_file = f"{code_dir}/meta/mm10.blacklist.bed"
#     output:
#         frag_file = expand("{o}/{{s}}_frag.h5ad", o = frag_dir),
#         frag_sum_file = expand("{o}/{{s}}_frag_sum.txt", o = frag_sum_dir),
#         qc_dlt_file = expand("{o}/{{s}}_qc_dlt.h5ad", o = qc_dlt_dir),
#         tag = touch(expand("{o}/{{s}}_qc_dlt.done", o = flag_dir))
#     log:
#         expand("{o}/{{s}}_qc_dlt.log", o = log_dir)
#     params:
#         sample = get_sample
#     threads: 1
#     resources:
#         walltime = "20:00:00",
#         queue = "hotel"
#     script:
#         f"{code_dir}/17.snapatac2/script/sa2.qc.dlt.py"

rule rm_dlt:
    input:
        qc_dlt_file = expand("{o}/{{s}}_qc_dlt.h5ad", o = qc_dlt_dir)
    output:
        snap_file = expand("{o}/{{s}}_rm_dlt.h5ad", o = rm_dlt_dir),
        tag = touch(expand("{o}/{{s}}_rm_dlt.done", o = flag_dir))
    log:
        expand("{o}/{{s}}_rm_dlt.log", o = log_dir)
    params:
        sample = get_sample
    threads: 1
    resources:
        walltime = "20:00:00",
        queue = "hotel"
    script:
        f"{code_dir}/17.snapatac2/script/sa2.rm.dlt.py"
    
rule merge:
    input:
        snap_files = expand("{o}/{s}_rm_dlt.h5ad",
                            o = rm_dlt_dir, s = samples)
    output:
        merge_snap = f"{out_dir}/merge_cemba_all.h5ad",
        tag = touch(f"{flag_dir}/merge_cemba_all.done")
    log:
        f"{log_dir}/merge_snap.log"
    threads: 4
    resources:
        walltime = "20:00:00",
        queue = "hotel"
    script:
        f"{code_dir}/17.snapatac2/script/sa2.merge.rmdlt.py"
    
