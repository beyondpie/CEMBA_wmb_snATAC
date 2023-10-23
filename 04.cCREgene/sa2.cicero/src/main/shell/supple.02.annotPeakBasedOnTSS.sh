projroot="/projects/ps-renlab2/szu/projects/CEMBA2"
allpeakFile="${projroot}/supple.07.peakcalling.allinone/mba.whole.sa2.final.peak.srt.bed"
tssBed="${projroot}/meta/gencode.vM23.gene.tssUpDn1k.bed"
outdir="${projroot}/04.cCREgene/sa2.cicero/src/main/resource"
allpeakOvlpTssBed1="${outdir}/mba.whole.sa2.peakOvlpTSS.bed"
allpeakOvlpTssBed2="${outdir}/mba.whole.sa2.peakOvlpTSS.proximal.distal.bed"
allpeakOvlpTssBed3="${outdir}/mba.whole.sa2.peakOvlpTSS.proximal.distal.ciceroPeakCoord.bed"


## annot distall, proximal peak
intersectBed -wao -f 0.5 -a ${allpeakFile} -b ${tssBed} > ${allpeakOvlpTssBed1}

awk 'BEGIN{FS=OFS="\t"}{if($5 != ".") {print $1,$2,$3,$4,"proximal",$11} else {print $1,$2,$3,$4,"distal","nan"}}' \
    ${allpeakOvlpTssBed1} \
    | sort -k1,1 -k2,2n | uniq > ${allpeakOvlpTssBed2}

awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$4,$5,$6}' ${allpeakOvlpTssBed2} \
    | sort -k1,1 > ${allpeakOvlpTssBed3}



