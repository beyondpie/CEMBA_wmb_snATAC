#!/bin/bash

sum_dir="/Users/szu/git-recipes/mouseBrainAtlas/CEMBA2/04.cCREgene/sa2.cicero/out/sa2pdcsum"
all_pdc_file="${sum_dir}/mba.whole.sa2subclass.merge.pdc.all"
distal_peak_nm="${sum_dir}/mba.whole.sa2subclass.all.distal.peaks.from.pdc.all.txt"
distal_peak_bed="${sum_dir}/mba.whole.sa2subclass.all.distal.peaks.from.pdc.all.bed"

awk 'BEGIN{FS=OFS="\t"}{NR>1}{print $5}' ${all_pdc_file} \
    | sort | uniq > ${distal_peak_nm}

awk -F '[:-]' 'BEGIN{OFS="\t"}{print $1,$2,$3}'  ${distal_peak_nm} \
    > ${distal_peak_bed}


