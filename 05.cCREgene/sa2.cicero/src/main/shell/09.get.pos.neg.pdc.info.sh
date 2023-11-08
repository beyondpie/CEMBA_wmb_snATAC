#!/bin/bash

proj_dir="/Users/szu/git-recipes/mouseBrainAtlas/CEMBA2"
cicero_dir="${proj_dir}/04.cCREgene/sa2.cicero/out/sa2pdcsum"
meta_col="sa2subclass"
allpeak="${proj_dir}/supple.07.peakcalling.allinone/mba.whole.sa2.final.peak.srt.bed"
pdc_bedpe=${cicero_dir}/mba.whole.${meta_col}.merge.bedpe.all
pdc_pair=${cicero_dir}/mba.whole.${meta_col}.merge.pdc.pair.all

# * function
sum_pos_or_neg_pdc () {
    local class=$1
    local cor_method=$2
    local class_pdc=${cicero_dir}/mba.whole.${meta_col}.${cor_method}.${class}.pdc.alignv1.tsv

    # outfiles
    local class_bedpe=${cicero_dir}/mba.whole.${meta_col}.${cor_method}.${class}.pdc.bedpe
    local class_pair=${cicero_dir}/mba.whole.${meta_col}.${cor_method}.${class}.pdc.pair
    local class_CREs=${cicero_dir}/mba.whole.${meta_col}.${cor_method}.${class}.pdc.CREs
    local class_bed=${cicero_dir}/mba.whole.${meta_col}.${cor_method}.${class}.pdc.CREs.bed
    
    # * process
    echo "$class under correlation method: ${cor_method}"
    
    # get bedpe
    # $8 is co-accessible score
    # the results have repeats since co-accessible are infered in subclass level and then
    # pool together
    join -1 1 -2 7 <(cut -f 1 ${class_pdc} | sort) <(sort -k7,7 ${pdc_bedpe})  -t$'\t' \
        | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1,$8,$9,$10}' \
        | sort -k1,1 -k2,2n | uniq > ${class_bedpe}
    # get pair
    join -1 1 -2 1 <(cut -f 1 ${class_pdc} | sort)  <(sort -k1,1 ${pdc_pair}) -t$'\t' \
        | sort -k1,1 -k2,2n | uniq > ${class_pair}
    # get CREs
    cut -f 1 ${class_pdc} | tr '|' '\t' | cut -f 2 | sort | uniq | sed '1d'> ${class_CREs}
    # get bed
    join -1 1 -2 4 ${class_CREs} <(sort -k4,4 ${allpeak}) -t$'\t' \
        | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1}' \
        | sort -k1,1 -k2,2n | uniq > ${class_bed}
    echo 'Done'
}

sum_pos_or_neg_pdc pos pearson
sum_pos_or_neg_pdc neg pearson


