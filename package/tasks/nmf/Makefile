.PHONY: nmf_subclassIntv2
nmf_subclassIntv2_encoder: nmf.Snakefile
	-mkdir -p $@
	cp config.yaml $@/config.yaml
	cp $< $@/Snakefile
	cp -R profile/. $@/profile
	cd $@ && \
  snakemake --config \
    system=encoder \
    tag=all.Intv2 \
    module=150 \
    n_rerun=2 \
    out=$@ \
    mat_pbyc_h5=data/cpm.cbyp.Intv2.h5 \
    peak_nm_file=data/peaks.Intv2.txt \
    cluster_nm_file=data/clusters.Intv2.txt \
    local_dir=package/tasks/nmf \
     -c 2 -p --snakefile Snakefile -R --rerun-incomplete

.PHONY: nmf_subclassIntv2_novlp_encoder
nmf_subclassIntv2_novlp_encoder: nmf.Snakefile
	-mkdir -p $@
	cp config.yaml $@/config.yaml
	cp $< $@/Snakefile
	cp -R profile/. $@/profile
	cd $@ && \
  snakemake --config \
    system=encoder \
    tag=all.Intv2 \
    mod_from=150 \
    mod_to=151 \
    mod_by=1 \
    n_rerun=2 \
    out=$@ \
    mat_pbyc_h5=data/cpm.cbyp.novlp.Intv2.h5 \
    peak_nm_file=data/peaks.novlp.Intv2.txt \
    cluster_nm_file=data/clusters.novlp.Intv2.txt \
    local_dir=package/tasks/nmf \
     -c 4 -p --snakefile Snakefile -R --rerun-incomplete


.PHONY: nmf_subclassIntv2_novlp_tscc
nmf_subclassIntv2_novlp_tscc: nmf.Snakefile
	-mkdir -p $@
	cp config.yaml $@/config.yaml
	cp $< $@/Snakefile
	cp -R profile/. $@/profile
	cd $@ && \
  snakemake --config \
    system=tscc \
    tag=novlp.Intv2 \
    module=150 \
    n_rerun=2 \
    out=$@ \
    mat_pbyc_h5=data/cpm.cbyp.novlp.Intv2.h5 \
    peak_nm_file=data/peaks.novlp.Intv2.txt \
    cluster_nm_file=data/clusters.novlp.Intv2.txt \
    local_dir=package/tasks/nmf \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete

.PHONY: nmf_subclassIntv2_all_tscc
nmf_subclassIntv2_all_tscc: nmf.Snakefile
	-mkdir -p $@
	cp config.yaml $@/config.yaml
	cp $< $@/Snakefile
	cp -R profile/. $@/profile
	cd $@ && \
  snakemake --config \
    system=tscc \
    tag=all.Intv2 \
    module=150 \
    n_rerun=2 \
    out=$@ \
    mat_pbyc_h5=data/cpm.cbyp.Intv2.h5 \
    peak_nm_file=data/peaks.Intv2.txt \
    cluster_nm_file=data/clusters.Intv2.txt \
    local_dir=package/tasks/nmf \
     -c 2 -p --snakefile Snakefile -R --rerun-incomplete --profile profile


clean:
	-rm -rf nmf_test
	-rm -rf nmf_test_local
	-rm -rf nmf_test_qsub
