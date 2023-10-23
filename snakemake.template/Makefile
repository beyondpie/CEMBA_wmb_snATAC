test_dir:
	-mkdir -p $@

.PHONY: test_snakemake
test_snakemake : test_dir Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    pmat_dir=$< \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile


submit_qsub:
	qsub pbs.demo.sh
