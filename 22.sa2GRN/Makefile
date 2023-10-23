.PHONY: celloracle

snakemake := /home/szu/miniforge3/envs/sa2/bin/snakemake
# snakemake := /home/szu/mambaforge/envs/sa22/bin/snakemake
# celloracle: src/main/pipeline/celloracle.Snakefile src/main/resource/config.yaml
# 	-mkdir -p build/$@
# 	cp $(word 2,$^) build/$@/config.yaml
# 	cp $< build/$@/Snakefile
# 	cd build/$@ && \
# 	${snakemake} --snakefile Snakefile  -R -c 30\
#      --config \
#      threads=1 \
#      debug=0 \
#      queue=glean \
#      walltime=8 \
#      system=mediator \
#      njobs=1 \
#     --rerun-triggers mtime \
#     --skip-script-cleanup \
#     --keep-going \
    # --profile pbs-torque-conda

# comment runGRN_baseGRNsc firstly in snakefile
celloracle: src/main/pipeline/celloracle.Snakefile src/main/resource/config.yaml
	-mkdir -p build/$@
	cp $(word 2,$^) build/$@/config.yaml
	cp $< build/$@/Snakefile
	cd build/$@ && \
	${snakemake} --snakefile Snakefile  -R -c 30\
     --config \
     threads=1 \
     debug=0 \
     queue=glean \
     walltime=8 \
     system=mediator \
     njobs=4 \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --keep-going \
    # --profile pbs-torque-conda
