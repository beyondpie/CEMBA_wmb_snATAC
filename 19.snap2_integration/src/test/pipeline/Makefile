.PHONY: test_R_snakemake test_intgn_snakemake

test_R_snakemake: R.Snakefile
	snakemake --snakefile $< -R -c 2\
    --rerun-triggers mtime \
    --skip-script-cleanup
	rm out/*.out

conda_path := /home/szu/mambaforge/envs/seurat/bin
test_intgn_snakemake: ../../main/pipeline/Seurat.Intgn.Snakefile
	-mkdir -p $@
	cp $< $@/Snakefile
	cd $@ && \
	snakemake --snakefile Snakefile -R -c 2 \
    --config \
    gp=nn \
    debug=1 \
    intgn_method='rpac,mnn' \
    allen_techs='10xv3' \
	  --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda

clean:
	-rm -rf out
	-rm -rf log
	-rm -rf .snakemake
