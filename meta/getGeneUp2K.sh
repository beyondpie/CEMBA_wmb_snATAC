#!/bin/bash

gtf="modified_gencode.vM23.primary_assembly.annotation.gtf"
bedfile="mouse.modified_gencode.vM23.bed"
geneUp2kfile="mmouse.modified.genecode.vM23.gene.up2k.bed"

awk 'BEGIN{FS=OFS="\t"}($3=="gene"){split($9,a,"\""); print $1,$4-1,$5,a[6],$6,$7}' ${gtf}\
    | sort -k1,1 -k2,2n > ${bedfile}

awk 'BEGIN{FS=OFS="\t"}{if($6=="+" && $2-2000>0){print $1,$2-2000,$3,$4,$5,$6}else if($6=="+" && $2-2000<0){print $1,0,$3,$4,$5,$6}else if($6=="-"){print $1,$2,$3+2000,$4,$5,$6}}' ${bedfile} > ${geneUp2kfile}
