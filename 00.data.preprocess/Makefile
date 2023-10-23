sa2_pp: snapatac2.qc.Snakefile
	snakemake -c 1 --config system=silencer \
    --snakefile $< -R --rerun-incomplete --profile profile


encoder_sa2_l1_embed: snapatac2.qc.Snakefile
	snakemake -c 4 --config system=silencer debug=1 \
    --snakefile $< --until snapatac2_l1_embed -R --rerun-incomplete

encoder_sa2_l1_knn_hora: snapatac2.qc.Snakefile
	snakemake -c 4 --config system=silencer debug=1 \
    knn_method=hora \
    --snakefile $< --until snapatac2_l1_knn -R --rerun-incomplete

encoder_sa2_l1_knn_exact: snapatac2.qc.Snakefile
	snakemake -c 4 --config system=silencer debug=1 \
    knn_method=exact \
    --snakefile $< --until snapatac2_l1_knn -R --rerun-incomplete

.PHONY: tscc_sa2_l1_embed
tscc_sa2_l1_embed: snapatac2.qc.Snakefile snapatac2.qc.config.yaml
	-mkdir -p $@
	cp $(word 2,$^) $@/$(word 2,$^)
	cp $< $@/$<
	cp -R profile/. $@/profile
	cd $@ && \
	snakemake -c 1 --config system=tscc debug=1 \
    --snakefile $< --until snapatac2_l1_embed -R \
    --rerun-incomplete --profile profile

.PHONY: tscc_sa2_l1_knn_exact
l1_knn_dir := tscc_sa2_l1_embed
tscc_sa2_l1_knn_exact: snapatac2.qc.Snakefile snapatac2.qc.config.yaml
	cp $(word 2,$^) ${l1_knn_dir}/$(word 2,$^)
	cp $< ${l1_knn_dir}/$<
	cp -R profile/. ${l1_knn_dir}/profile
	cd ${l1_knn_dir} && \
	snakemake -c 1 --config \
    system=tscc debug=1 \
    knn_method=exact \
    --snakefile $< --until snapatac2_l1_knn -R \
    --rerun-incomplete --profile profile


# test
test_sa2_pp: snapatac2.qc.Snakefile
	snakemake -c 1 --config system=imac debug=1\
    --snakefile $< -R --rerun-incomplete

test_sa2_pp_tscc: snapatac2.qc.Snakefile
	snakemake -c 1 --config system=tscc debug=1 \
    --snakefile $< -R --rerun-incomplete

test_sa2_l1_embed: snapatac2.qc.Snakefile
	snakemake -c 2 --config system=imac debug=1 \
    --snakefile $< --until snapatac2_l1_embed -R --rerun-incomplete

test_sa2_l1_knn: snapatac2.qc.Snakefile
	snakemake -c 2 --config system=imac debug=1 \
    --snakefile $< --until snapatac2_l1_knn -R --rerun-incomplete

.PHONY: clean
clean:
	-rm snapatac2_l1_embed..03*
	-rm snapatac2_pp.sample=*.o3*
	-rm -rf .tmp*
