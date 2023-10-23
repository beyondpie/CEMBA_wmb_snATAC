proj_dir := /projects/ps-renlab2/szu/projects/CEMBA2
snakemake := /home/szu/mambaforge/envs/sa2/bin/snakemake
oasis_dir := /oasis/tscc/scratch/szu/projects/CEMBA2

# * get bedfiles
.PHONY: neu_all_bed neu_rep_bed

neu_all_bed: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 64 \
	   --config \
      gp=neuron \
      sys=tscc \
      queue=home \
      walltime=24 \
      thread=5\
      debug=0 \
      bedtag='all' \
      nmin_macs2=400 \
      macs2=/home/szu/mambaforge/bin/macs2 \
      cluster2size_fnm=neuron_L4pc2sizes_cca-k50-cl_v1.csv \
    --until getBedFile \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda 

neu_rep_bed: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 64 \
	   --config \
      gp=neuron \
      sys=tscc \
      queue=home \
      walltime=24 \
      thread=5\
      debug=0 \
      bedtag='biorep,pseudorep' \
      nmin_macs2=400 \
      macs2=/home/szu/mambaforge/bin/macs2 \
      out_macs2=/oasis/tscc/scratch/szu/projects/CEMBA2/18.snap2_peakcalling \
      cluster2size_fnm=neuron_L4pc2sizes_cca-k50-cl_v1.csv \
    --until getBedFile \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda

.PHONY: encoder_fix_neuron_bed encoder_nn_bed
encoder_fix_neuron_bed: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 3 \
	   --config \
      gp=neuron \
      sys=encoder \
      queue=home \
      walltime=24 \
      thread=1\
      debug=1 \
      bedtag='all' \
      nmin_macs2=200 \
      macs2=/home/szu/mambaforge/envs/sa2/bin/macs2 \
      out_macs2=/projects/ps-renlab2/szu/projects/CEMBA2/18.snap2_peakcalling \
      cluster2size_fnm=neuron_L4pc2sizes_cca-k50-cl_v1.csv \
    --until getBedFile \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup

encoder_nn_bed: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 3 \
	   --config \
      gp=nn \
      sys=encoder \
      queue=home \
      walltime=24 \
      thread=1\
      debug=0 \
      bedtag='all,biorep,pseudorep' \
      nmin_macs2=200 \
      macs2=/home/szu/mambaforge/envs/sa2/bin/macs2 \
      out_macs2=/projects/ps-renlab2/szu/projects/CEMBA2/18.snap2_peakcalling \
      cluster2size_fnm=nn_L4pc2sizes_cca-k49-cl_v1.csv \
    --until getBedFile \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup


# * peakcallinng
.PHONY: tscc_neu_macs2 tscc_nn_macs2
tscc_neu_macs2: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 30 \
	   --config \
      gp=neuron \
      sys=tscc \
      queue=glean \
      walltime=1 \
      thread=5\
      debug=0 \
      bedtag='all,biorep,pseudorep' \
      nmin_macs2=200 \
      macs2=/home/szu/mambaforge/envs/sa22/bin/macs2 \
      out_macs2=/oasis/tscc/scratch/szu/projects/CEMBA2/18.snap2_peakcalling \
      cluster2size_fnm=neuron_L4pc2sizes_cca-k50-cl_v1.csv \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda 

tscc_nn_macs2: src/main/pipeline/snap2.peakcalling.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 30 \
	   --config \
      gp=nn \
      sys=tscc \
      queue=condo \
      walltime=3 \
      thread=5\
      debug=0 \
      bedtag='all,biorep,pseudorep' \
      nmin_macs2=200 \
      macs2=/home/szu/mambaforge/envs/sa22/bin/macs2 \
      out_macs2=/oasis/tscc/scratch/szu/projects/CEMBA2/18.snap2_peakcalling \
      cluster2size_fnm=nn_L4pc2sizes_cca-k49-cl_v1.csv \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda

# * merge peak
.PHONY: merge_peak
merge_peak: src/main/pipeline/mergePeak.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 20 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup
    # --profile pbs-torque-conda

rsc_dir := src/main/resource
.PHONY: all_pL4s
all_pL4s: ${rsc_dir}/neuron_L4pc2sizes_cca-k50-cl_v1.csv ${rsc_dir}/nn_L4pc2sizes_cca-k49-cl_v1.csv
	-rm ${rsc_dir}/all_pl4s.txt
	awk -F ',' '{print "neuron."$$1}' $(word 1,$^) > all_pL4s.txt && \
	awk -F ',' '{print "nn."$$1}' $(word 2,$^) >> all_pL4s.txt && \
	cp all_pL4s.txt > ${rsc_dir}/all_pL4s.txt

# intersect_mergepeak: src/main/shell/intersect_mergepeak.sh
# 	parallel -j 8 bash $< out/tscc/mba.whole.union.peak.bed  \
#     out/tscc/itermerge/{}.filterNfixed.peakset \
#     out/tscc/unionpeak/{}.union.peak.bed :::: ${rsc_dir}/all_pL4s.txt

# * get pmat
.PHONY: get_pmat
get_pmat : src/main/pipeline/getsa2pmat.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 9 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    #--profile pbs-torque-conda

# # * fitscmodel of background for frac
# .PHONY: fitbg
fitbg : src/main/pipeline/scfilter.Snakefile
	-mkdir -p build/$@
	cp $< build/$@/Snakefile
	cd build/$@ && ${snakemake} --snakefile Snakefile -R -c 9 \
    --keep-incomplete \
    --rerun-triggers mtime \
    --skip-script-cleanup \
    --profile pbs-torque-conda
