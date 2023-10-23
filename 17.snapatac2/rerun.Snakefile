import os

configfile: "config.yaml"

system = config["system"]
proj_dir = config["project_dir"][system]
work_dir = config["work_dir"][system]
out_dir = f'{work_dir}/{config["out_dir"]}'
log_dir = f"{out_dir}/log"
flag_dir = f"{out_dir}/flag"
os.makedirs(out_dir, exist_ok = True)
os.makedirs(log_dir, exist_ok = True)
os.makedirs(flag_dir, exist_ok = True)

conda_env = config["conda"][system]

cll = config["clustering_level"]
print("Clustering on level: ", cll)

# NOTE: remote clustering level info temporally in pcm_file loading
pcm_file = f'{proj_dir}/{config["pre_clustering_meta"]}'
with open(pcm_file, 'r') as f:
    lines = [l.strip() for l in f.readlines()]
pre_id2size = {l.split(",")[0]: int(l.split(",")[1]) for l in lines}

# clustering: queue glean with ncore: 4
# pre_clusters = pre_id2size.keys()
pre_clusters = [1,2,3,11,19,32,34,41,46,52,53,56]
max_united_size = config['max_united_size']
print(f"Max size for united clustering: {max_united_size}")


blacklist_file = f'{proj_dir}/{config["blacklist_file"]}'
barcode2id_file = f'{proj_dir}/{config["barcode2id_file"]}'
sample2fragment_file = f'{proj_dir}/{config["sample2fragment_file"]}'

for f in [blacklist_file, barcode2id_file, sample2fragment_file]:
    if not os.path.exists(f):
        raise RuntimeError(f"{f} is not found.")

def get_clusterid(wildcards):
    return wildcards.i

def get_ppn(wildcards, attempt):
    return attempt * 8

def get_walltime_embed(wildcards, attempt):
    wt = attempt * 30
    return f"{wt}:00:00"

def get_walltime_knn(wildcards, attempt):
    wt = attempt * 30
    return f"{wt}:00:00"

def get_walltime_leiden(wildcards, attempt):
    wt = attempt * 30
    return f"{wt}:00:00"

def get_walltime_unite(wildcards, attempt):
    wt = attempt * 24
    return f"{wt}:00:00"

# all move to hotel
def get_queue(wildcards, attempt):
    wt = attempt * 4
    if wt > 8:
        return "hotel"
    else:
        return "hotel"
    
def get_clustering_input(wildcards):
    if pre_id2size[wildcards.i] > max_united_size:
        print("Run rules in seperated embed, knn way.")
        return f"{flag_dir}/leiden_of_{cll}_{wildcards.i}.done"
    else:
        print("Run rules in united way.")
        return f"{flag_dir}/united_of_{cll}_{wildcards.i}.done"
    
embed_params = config["embed"]
knn_params = config["knn"]
leiden_params = config["leiden"]
umap_params = config["umap"]

rule all:
    input:
        expand("{d}/clustering_of_{l}_{i}.done",
               d = flag_dir, l = cll, i = pre_clusters)

rule anndataset:
    input:
        sample2fragment = sample2fragment_file
    output:
        f"{out_dir}/cemba_all_anndataset.hdf5"
    log:
        f"{log_dir}/cemba_all_anndataset.log"
    threads:1
    resources:
        walltime = "05:00:00",
        queue = "hotel"
    script:
        "script/sa2.pre.anndataset.py"
        
rule embed:
    input:
        blacklist = blacklist_file,
        barcode2id = barcode2id_file,
        cemba_anndataset = rules.anndataset.output[0]
    output:
        snap_file = expand("{o}/{nm}_{{i}}_mult.h5ad",
               o = out_dir, nm = embed_params["name"]),
        tag = touch(expand("{d}/embed_of_{l}_{{i}}.done",
                           d = flag_dir, l = cll))
    log:
        expand("{d}/embed_of_{l}_{{i}}.log",
               d = log_dir, l = cll)
    threads: get_ppn
    params:
        embed = embed_params,
        cluster_id = get_clusterid,
        clevel = cll
    retries: 2
    resources:
        walltime = get_walltime_embed,
        queue = get_queue
    script:
        "script/sa2.embed.py"

rule knn:
    input:
        snap_file = expand("{o}/{nm}_{{i}}_mult.h5ad",
                           o = out_dir, nm = embed_params["name"]),
        tag = expand("{d}/embed_of_{l}_{{i}}.done", d = flag_dir, l = cll)
    output:
        touch(expand("{d}/knn_of_{l}_{{i}}.done",
                     d = flag_dir, l = cll))
    log:
        expand("{d}/knn_of_{l}_{{i}}.log",
               d = log_dir, l = cll)
    retries: 2
    params:
        knn = knn_params
    threads: get_ppn
    resources:
        walltime = get_walltime_knn,
        queue = get_queue
    script:
        "script/sa2.knn.py"

rule leiden:
    input:
        tag= expand("{d}/knn_of_{l}_{{i}}.done", d = flag_dir, l = cll),
        snap_file = expand("{o}/{nm}_{{i}}_mult.h5ad",
                           o = out_dir, nm = embed_params["name"])
    output:
        touch(expand("{d}/leiden_of_{l}_{{i}}.done",
                     d = flag_dir, l = cll))
    log:
        expand("{d}/leiden_of_{l}_{{i}}.log",
               d = log_dir, l = cll)
    retries: 2
    params:
        leiden = leiden_params,
        umap = umap_params,
        knn = knn_params,
        embed = embed_params,
        cluster_id = get_clusterid,
        clevel = cll
    threads: get_ppn
    resources:
        walltime = get_walltime_leiden,
        queue = get_queue
    script:
        "script/sa2.leiden.py"

# TODO: add summary
# - draw umap and clusters with different resolution
# - draw sample balance
rule united:
    input:
        blacklist = blacklist_file,
        barcode2id = barcode2id_file,
        cemba_anndataset = rules.anndataset.output[0]
    output:
        snap_file = expand("{o}/{nm}_{{i}}_unite.h5ad",
               o = out_dir, nm = embed_params["name"]),
        flag = touch(expand("{d}/united_of_{l}_{{i}}.done",
                            d = flag_dir, l = cll))
    log:
        expand("{d}/united_of_{l}_{{i}}.log", d = log_dir, l = cll)
    retries: 3
    params:
        embed = embed_params,
        knn = knn_params,
        leiden = leiden_params,
        umap = umap_params,
        cluster_id = get_clusterid,
        clevel = cll
    threads: get_ppn
    resources:
        walltime = get_walltime_unite,
        # May need queue as hotel for united.
        queue = "hotel"
    script:
        "script/sa2.united.py"

# checkpoint
rule clustering:
    input:
        flag = get_clustering_input
    output:
        touch(expand("{d}/clustering_of_{l}_{{i}}.done", d = flag_dir, l = cll))
    threads: get_ppn
    resources:
        walltime = "01:00:00",
        queue = "glean"

        
