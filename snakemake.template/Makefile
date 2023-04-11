pmatTest:
	-mkdir -p $@

.PHONY: test_pmat
test_pmat : pmatTest cellpmat.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    group_file=04.cCREgene/test.group.list \
    pmat_dir=$< \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

pmatAllen:
	-mkdir -p $@
.PHONY: allen_pmat
allen_pmat : pmatAllen cellpmat.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    pmat_dir=$< \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

pmatL3:
	-mkdir -p $@
.PHONY: L3_pmat
L3_pmat : pmatL3 cellpmat.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    pmat_dir=$< \
    group_file=meta/L3.only.csv \
    meta_col=L3 \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile


ciceroTest:
	-mkdir -p $@
.PHONY: test_cicero
test_cicero: ciceroTest cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    group_file=04.cCREgene/test.group.list \
    cicero_dir=$< \
    debug=1 \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

ciceroAllen:
	-mkdir -p $@
.PHONY: runCicero
runCicero: ciceroAllen cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
    cicero_dir=$< \
	pmat_dir=pmatAllen \
    debug=0 \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

.PHONY: testFilterByShuf
testFilterByShuf: ciceroTest cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=test_cicero \
     group_file=04.cCREgene/test.group.list \
     debug=1 \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete


.PHONY: runFilterByShufAllen
runFilterByShufAllen: ciceroAllen cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=ciceroAllen \
     debug=0 \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

.PHONY: testFilterDistalProximalConns
testFilterDistalProximalConns: ciceroTest cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=test_cicero \
     group_file=04.cCREgene/test.group.list \
     -c 2 -p --snakefile Snakefile -R --rerun-incomplete

.PHONY: filterDistalProximalConns
filterDistalProximalConns: ciceroAllen cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=ciceroAllen \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile

.PHONY: testMergeDistalProximalConns
testMergeDistalProximalConns: ciceroTest cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=test_cicero \
     -c 3 -p --snakefile Snakefile -R --rerun-incomplete

.PHONY: mergeDistalProximalConns
mergeDistalProximalConns: ciceroAllen cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=ciceroAllen \
     -c 3 -p --snakefile Snakefile -R --rerun-incomplete

.PHONY: getCor
getCor: ciceroAllen cicero.Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config \
     cicero_dir=ciceroAllen \
      -c 3 -p --snakefile Snakefile -R --rerun-incomplete
