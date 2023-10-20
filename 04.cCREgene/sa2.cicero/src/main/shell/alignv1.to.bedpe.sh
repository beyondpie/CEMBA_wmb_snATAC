#!/bin/bash

all_pdc_bedpe=$1
pdc_alignv1=$2

join -1 1 -2 7 <(cut -f 1 ${pdc_alignv1} | sort) <(sort -k7,7 ${all_pdc_bedpe})  -t$'\t' \
    | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1,$8,$9,$10}' \
    | sort -u -k1,1 -k2,2n -k8,8nr 
