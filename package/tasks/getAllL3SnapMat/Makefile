mat := vM16gmat

define screenSnakemake
	screen -dmS ${1} snakemake -c 1 -p --snakefile Snakefile --profile profile -R --rerun-incomplete
endef

run: Snakefile.template L2GroupAll.csv
	cp $< ${mat}/Snakefile
	cp -R profile.template/. ${mat}/profile
	cp $(word 2,$^) ${mat}/group.csv
	cd ${mat} && $(call screenSnakemake,${mat})
