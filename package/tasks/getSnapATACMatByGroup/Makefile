mat := L2vM16gmat

define screenSnakemake
	screen -dmS ${1} snakemake -c 1 -p --snakefile Snakefile --profile profile -R --rerun-incomplete
endef

run: Snakefile.template L2GroupAll.csv config.json.template
	-mkdir -p ${mat}
	cp $< ${mat}/Snakefile
	cp -R profile.template/. ${mat}/profile
	cp L2GroupAll.csv ${mat}/group.csv
	cp config.json.template ${mat}/config.json
	cd ${mat} && $(call screenSnakemake,${mat})

test: Snakefile.template L2GroupTest.csv config.json.template
	-mkdir -p ${mat}
	cp $< ${mat}/Snakefile
	cp -R profile.template/. ${mat}/profile
	cp L2GroupTest.csv ${mat}/group.csv
	cp config.json.template ${mat}/config.json
	cd ${mat} && $(call screenSnakemake,${mat})

l1mat := L1vM16gmat
runL1: Snakefile.template L1Group.csv configL1.json.template
	-mkdir -p ${l1mat}
	cp $< ${l1mat}/Snakefile
	cp -R profile.template/. ${l1mat}/profile
	cp $(word 2,$^) ${l1mat}/group.csv
	cp $(word 3,$^) ${l1mat}/config.json
	cd ${l1mat} && $(call screenSnakemake,${l1mat})

multiL2 := L2MultvM16gmat
runMultiL2: Snakefile.template L2MultiGroup.csv configMultiGroup.json.template
	-mkdir -p ${multiL2}
	cp $< ${multiL2}/Snakefile
	cp -R profile.template/. ${multiL2}/profile
	cp $(word 2,$^) ${multiL2}/group.csv
	cp $(word 3,$^) ${multiL2}/config.json
	cd ${multiL2} && $(call screenSnakemake,${multiL2})
