snakemake_encoder := /home/szu/mambaforge/envs/sa2/bin/snakemake

# * quality control and doublet removal
.PHONY: qc_dlt filter_qc_dlt
qc_dlt: sa2.qc.dlt.Snakefile
	-mkdir -p $@
	cp $< $@/$<
	cd $@ && \
	${snakemake_encoder} --cores 8 \
    --snakefile $< -R \
    --rerun-triggers mtime
    # --profile pbs-torque-conda \

python := /home/szu/mambaforge/envs/sa2/bin/python
filter_qc_dlt: sa2.bmat.dlt.py
	parallel --bar -k -a ../meta/mba.whole.sample.lst ${python} $<

# * gmat
# need to activate conda env firstly
.PHONY: gmat
gmat: sa2.gmat.Snakefile
	snakemake --cores 8 \
    --snakefile $< -R \
    --rerun-triggers mtime

# * L1-level clustering
.PHONY: L1_clustering_umap L1_encoder
# run on tscc
L1_clustering_umap: script/sa2.clustering.umap.py
	python $< L1_encoder L0 L1_encoder/umap 30000 0

L1_encoder: Snakefile
	${snakemake_encoder} --cores 8 --snakefile $< -R \
     --config \
     system=encoder \
	   out_dir=$@ \
	   clustering_level=L0 \
     pre_clustering_meta=17.snapatac2/resource/sa2_L0_cluster2size.csv \
		 sample2fragment_file=17.snapatac2/resource/bmatfile_L0_condo.csv \
		 barcode2id_file=17.snapatac2/resource/barcode2id_L0.csv \
     embed_nfeat=5000000 \
     embed_ncomp=50 \
     embed_name=nfeat-all_nsample-all_nc50 \
		 max_united_size=60000\
     --rerun-triggers mtime

L1_add_barcode_unite: script/supple.sa2.add.barcode.to.unite.clustering.py
	python $<

# * L2-level clustering
snakemake_tscc := /home/szu/mambaforge/envs/sa2/bin/snakemake
.PHONY: L2_dlt2_encoder L2_dlt2_UMAP
# will have 300 Mem sometimes
# In total, takes about 12 hours for 37 clusters
L2_dlt2_encoder: Snakefile
	${snakemake_tscc} --cores 8 \
    --snakefile $< -R \
    --config \
    system=encoder \
    out_dir=$@ \
    clustering_level=L1 \
		pre_clustering_meta=17.snapatac2/resource/sa2_dlt2_L1_cluster2size.csv \
		barcode2id_file=17.snapatac2/resource/sa2_dlt2_L1_barcode2id.csv \
		max_united_size=100000 \
    embed_nfeat=500000 \
    embed_ncomp=50\
    embed_name=nfeat-top_nc50 \
    --rerun-triggers mtime

# L2 UMAP
L2_dlt2_UMAP: script/sa2.clustering.umap.py
	parallel python $< result/L2_sa2_dlt2_sum L1 result/L2_UMAP 30000  ::: {0..36}

# * L3 clustering
.PHONY: L3_dlt2_encoder L3_dlt2_UMAP
L3_dlt2_encoder: Snakefile
	${snakemake_tscc} --cores 8 \
    --snakefile $< -R \
    --config \
    system=encoder \
    out_dir=$@ \
    clustering_level=L2 \
		pre_clustering_meta=17.snapatac2/resource/sa2_dlt2_L2_cluster2size.csv \
		barcode2id_file=17.snapatac2/resource/sa2_dlt2_L2_barcode2id.csv \
		max_united_size=100000 \
    embed_nfeat=500000 \
    embed_ncomp=30\
    embed_name=nfeat-top_nc30 \
    --rerun-triggers mtime

L3_dlt2_cids.txt : result/L3_sa2_dlt2_pkl
	ls -lA $< | awk '{print $$9}' | awk -F "[_.]" '{print $$4}' | sed '/^$$/d'> $@

python := ~/mambaforge/envs/sa2/bin/python
L3_dlt2_UMAP: script/sa2.clustering.umap.py L3_dlt2_cids.txt
	cat $(word 2,$^) | parallel ${python} $< result/L3_sa2_dlt2_pkl L2 result/L3_sa2_dlt2_UMAP 30000

# * L4 clustering
.PHONY: L4_dlt2_encoder L4_dlt2_UMAP
L4_dlt2_encoder : Snakefile
	${snakemake_tscc} --cores 4 \
		--snakefile $< -R \
		--config \
		system=encoder \
		out_dir=$@ \
		clustering_level=L3 \
		pre_clustering_meta=17.snapatac2/resource/sa2_dlt2_L3forL4_cluster2size.csv \
		barcode2id_file=17.snapatac2/resource/sa2_dlt2_L3forL4_barcode2id.csv \
		max_united_size=200000 \
		embed_nfeat=500000 \
		embed_ncomp=30\
		embed_name=nfeat-top_nc30 \
		--rerun-triggers mtime

L4_dlt2_cids.txt : result/L4_sa2_dlt2_pkl
	ls -lA $< | awk '{print $$9}' | awk -F "[_.]" '{print $$4}' | sed '/^$$/d'> $@

python := ~/mambaforge/envs/sa2/bin/python
L4_dlt2_UMAP: script/sa2.clustering.umap.py L4_dlt2_cids.txt
	cat $(word 2,$^) | parallel ${python} $< result/L4_sa2_dlt2_pkl L3 result/L4_sa2_dlt2_UMAP 30000


clean:
	-rm test*.log
	-rm -rf test_qc_dlt
	-rm -rf test_tscc_unite
	-rm -rf .tmp*
	-rm L3_dlt2_cids.txt
