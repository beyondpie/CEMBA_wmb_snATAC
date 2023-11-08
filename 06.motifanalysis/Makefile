.PHONY: ppdcpeak_motif
snakemake := /home/szu/miniforge3/bin/snakemake
# * run non-negative matrix factorization 
ppdcpeak_motif: src/main/pipeline/motif.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && \
  ${snakemake} --config \
    system=mediator \
    tag=ppdcpeak \
    out=ppdcpeak \
    njob=4 \
    module=54 \
     -c 20 \
     --snakefile Snakefile -R --rerun-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup  \
     # --profile pbs-torque-conda

clean_homer_tag:
	-rm out/ppdcpeak/flag/nmf.ppdcpeak.homer.*.n*.done
