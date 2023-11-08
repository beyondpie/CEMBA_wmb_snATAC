.PHONY: install_seuratv5

# * neuron + 10xv3 integration
.PHONY: neuron_10xv3 neuron_10xv3_encoder
snakemake := /home/szu/mambaforge/envs/seurat/bin/snakemake
neuron_10xv3: src/main/pipeline/Seurat.Intgn.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
	${snakemake} --snakefile Snakefile -R -c 32 \
    --config \
    gp=neuron \
    debug=0 \
    intgn_method='rpca' \
    allen_techs='10xv3' \
    feature='TF,markers,k8' \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda

# * neuron 10xv3 transfer data
## encoder
.PHONY: neu_tf_ds50_encoder
neu_tf_ds50_encoder: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 1 \
     --config \
     gp=neuron \
     sys=encoder \
     queue=encoder\
     debug=0 \
     intgn_method='cca' \
     allen_techs='10xv3' \
     feature='k8' \
     k_anchors='50,70,100' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=150\
     atac_ds_max=151\
     allen_ds_min=50\
     allen_ds_max=51\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup 

## tscc
.PHONY: neu_tf_ds50_tscc
neu_tf_ds50_tscc: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 64 \
     --config \
     gp=neuron \
     sys=tscc \
     queue=home \
     debug=0 \
     intgn_method='cca' \
     allen_techs='allen' \
     feature='k8' \
     k_anchors='100,150,200' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=150\
     atac_ds_max=151\
     allen_ds_min=50\
     allen_ds_max=51\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup \
     --profile pbs-torque-conda

.PHONY: imnneu_tf_ds50_tscc
# only neuronal cells from both data in Allen IMN top regions 
# no IMN from nn involved 
imnneu_tf_ds50_tscc: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 50 \
     --config \
     gp=imneuron \
     sys=tscc \
     queue=glean \
     debug=0 \
     intgn_method='cca' \
     allen_techs='allen' \
     feature='k8' \
     k_anchors='50' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=150\
     atac_ds_max=151\
     allen_ds_min=50\
     allen_ds_max=51\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup
     # --profile pbs-torque-conda

.PHONY: imnneu2_tf_ds50_tscc
# due to the bug of Seuratv5 merge
# https://github.com/satijalab/seurat/issues/7191
# only data.1 are used, so equals to the one above
imnneu2_tf_ds50_tscc: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 50 \
     --config \
     gp=imneuron \
     sys=tscc \
     queue=glean \
     debug=0 \
     intgn_method='cca' \
     allen_techs='nneuronallen' \
     feature='k8' \
     k_anchors='50' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=150\
     atac_ds_max=151\
     allen_ds_min=50\
     allen_ds_max=51\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup
     # --profile pbs-torque-conda

.PHONY: imnneu3_tf_ds50_tscc
# due to memory issue, directly downsample data
# and both IMN and neuronal cells are included from allen
# in total, about 200K cells are used
# needs > 167G RAM in CCA
# about 2 - 3 hours to finish on a 61 cores node
# it sometimes use 10 CPUs or more sometimes during CCA
# but usually 1 CPUs
imnneu3_tf_ds50_tscc: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 50 \
     --config \
     gp=imneuron \
     sys=tscc \
     queue=glean \
     debug=0 \
     intgn_method='cca' \
     allen_techs='dp50allen' \
     feature='k8' \
     k_anchors='50' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=150\
     atac_ds_max=151\
     allen_ds_min=50\
     allen_ds_max=51\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup
     # --profile pbs-torque-conda

.PHONY: nn_tf_ds200_tscc
nn_tf_ds500_tscc: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 20 \
     --config \
     gp=nn \
     sys=tscc \
     queue=glean \
     debug=0 \
     intgn_method='cca' \
     allen_techs='allen' \
     feature='k8' \
     k_anchors='49,100,150' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=300\
     atac_ds_max=301\
     allen_ds_min=500\
     allen_ds_max=501\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup \
     --profile pbs-torque-conda


.PHONY: dotplot_neu_tf dotplot_nn_tf
#dotplot_neu_tf: src/main/R/analyzetf.R
#	parallel -j 5 Rscript --vanilla $< \
#  ::: 5 30 ::: vf markers k8 TF ::: 10xv3 allen ::: rpca cca

dotplot_neu_tf: src/main/R/analyzetf.R
	parallel -j 2 Rscript --vanilla $< ::: 50 70 ::: k8 ::: allen ::: cca ::: neuron

dotplot_nn_tf: src/main/R/analyzetf.R
	parallel -j 2 Rscript --vanilla $< ::: 49 100 150 ::: k8 ::: allen ::: cca ::: nn

dotplot_imn_tf: src/main/R/analyzetf.R
	Rscript --vanilla $< 50 k8 dp50allen cca imneuron

lightsum_neu_tf: src/main/R/lightSummaryOfTransferLabel.R
	parallel -j 2 Rscript --vanilla $< ::: 30 50 ::: neuron ::: allen

lightsum_nn_tf: src/main/R/lightSummaryOfTransferLabel.R
	parallel -j 2 Rscript --vanilla $< ::: 49 50 ::: nn ::: allen

lightsum_imn_tf: src/main/R/lightSummaryOfTransferLabel.R
	Rscript --vanilla $< 50 imneuron dp50allen


# * test
.PHONY: test_tf
test_tf: src/main/pipeline/Seurat.TransferLabel.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
   ${snakemake} --snakefile Snakefile -R -c 4 \
     --config \
     gp=nn \
     sys=tscc \
     queue=glean \
     debug=0 \
     intgn_method='cca' \
     allen_techs='10xv3' \
     feature='vf' \
     allen_as_ref=1\
     ds_seu=1\
     atac_ds_min=20\
     atac_ds_max=21\
     allen_ds_min=10\
     allen_ds_max=11\
     --keep-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup \
     --profile pbs-torque-conda

# * clean up
clean_script:
	-rm build/*/.snakemake/scripts/*

clean:
	-rm -rf build/test_tf
