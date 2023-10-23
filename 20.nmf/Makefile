.PHONY: allpeak_nmf novlpeak_nmf ppdcpeak_nmf
# * run non-negative matrix factorization 

allpeak_nmf: src/main/pipeline/nmf.Snakefile src/main/resource/config.yaml
	-mkdir -p build/$@
	cp $(word 2,$^) build/$@/config.yaml
	cp $< build/$@/Snakefile
	cd build/$@ && \
  snakemake --config \
    system=encoder \
    tag=allpeak \
    out=out/allpeak_nmf \
    mod_from=150 mod_to=152 mod_by=2 n_rerun=3 \
    subclass_order_meta=data/sa2.subclass.srt.txt \
    peak_nm_file=data/sa2.final.peak.nms.txt \
    cluster_nm_file=data/sa2.allpeak.subclass.nms.txt \
    mat_pbyc_h5=data/sa2.allpeak.cbyp.mat.h5 \
     -c 10 \
     --snakefile Snakefile -R --rerun-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup  \
     # --profile pbs-torque-conda


novlpeak_nmf: src/main/pipeline/nmf.Snakefile src/main/resource/config.yaml
	-mkdir -p build/$@
	cp $(word 2,$^) build/$@/config.yaml
	cp $< build/$@/Snakefile
	cd build/$@ && \
  snakemake --config \
    system=encoder \
    tag=novlpeak \
    out=out/novlppeak_nmf \
    mod_from=150 mod_to=152 mod_by=2 n_rerun=3 \
    subclass_order_meta=data/sa2.subclass.srt.txt \
    peak_nm_file=data/sa2.novlpDHS.peak.nms.txt \
    cluster_nm_file=data/sa2.novlpDHS.subclass.nms.txt \
    mat_pbyc_h5=data/sa2.novlpDHS.cbyp.mat.h5 \
     -c 10 \
     --snakefile Snakefile -R --rerun-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup  \
     # --profile pbs-torque-conda



ppdcpeak_nmf: src/main/pipeline/nmf.Snakefile src/main/resource/config.yaml
	-mkdir -p build/$@
	cp $(word 2,$^) build/$@/config.yaml
	cp $< build/$@/Snakefile
	cd build/$@ && \
  snakemake --config \
    system=encoder \
    tag=ppdcpeak \
    out=out/ppdcpeak_nmf \
    mod_from=54 mod_to=56 mod_by=2 n_rerun=3 \
    subclass_order_meta=data/sa2.subclass.srt.txt \
    peak_nm_file=data/sa2.ppdc.peak.nms.txt \
    cluster_nm_file=data/sa2.ppdc.subclass.nms.txt \
    mat_pbyc_h5=data/sa2.ppdc.cbyp.mat.h5 \
     -c 10 \
     --snakefile Snakefile -R --rerun-incomplete \
     --rerun-triggers mtime \
     --skip-script-cleanup  \
     # --profile pbs-torque-conda
