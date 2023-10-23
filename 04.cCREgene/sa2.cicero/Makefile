.PHONY: run_cicero pdc_cicero

snakemake := /home/szu/miniforge3/bin/snakemake
# snakemake := /home/szu/mambaforge/envs/seurat/bin/snakemake
run_cicero: src/main/pipeline/runCicero.Snakefile
	${snakemake} --snakefile $< -R -c 9 \
    --config \
    threads=8 \
    debug=0 \
    queue=hotel \
    walltime=16 \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda

pdc_cicero: src/main/pipeline/pdc.Snakefile src/main/resource/config.yaml
	-mkdir -p build/$@
	cp $(word 2,$^) build/$@/config.yaml
	cp $< build/$@/Snakefile
	cd build/$@ && \
  ${snakemake} -c 20 --config \
    threads=1 \
    debug=0 \
    queue=glean \
    walltime=2 \
		--snakefile Snakefile -R --rerun-incomplete \
		--rerun-triggers mtime \
		--skip-script-cleanup  \
    # --profile pbs-torque-conda

getppdc_bedpe: src/main/shell/09.get.pos.neg.pdc.info.sh
	bash $<
