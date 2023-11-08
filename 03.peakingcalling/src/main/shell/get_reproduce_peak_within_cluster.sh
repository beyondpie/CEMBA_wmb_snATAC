#!/bin/bash
# example:
# bash 02d.naiveoverlap.sh -i c1 \
#  -p /projects/ps-renlab2/kaw033/CEMBAv2/CEMBA_data/snapatac2 \
#  -o /projects/ps-renlab2/kaw033/CEMBAv2/CEMBA_data/snapatac2 

while getopts ':i:p:o:h' opt; do
  case "$opt" in
    i)
      cl="$OPTARG"
      echo "Processing option 'i' with '${OPTARG}' argument"
      ;;
    p)
      from_dir="$OPTARG"
      echo "Processing option 'p' with '${OPTARG}' argument"
      ;;
    o)
      out_dir="$OPTARG"
      echo "Processing option 'o' with '${OPTARG}' argument"
      ;;
    h)
      echo "Usage: $(basename $0) [-i cl] [-p from_dir] [-o out_dir]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-i cl] [-p from_dir] [-o out_dir]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-i cl] [-p from_dir] [-o out_dir]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


if [ ! -d "${out_dir}" ]; then
  mkdir -p ${out_dir}
fi

#BLACKLIST=/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz

function reuse_pooled_peak ()
{
    local r=$1
    if [ ! -f $1 ]; then
        r=$2
    fi
    echo "$r"
}

POOLED_PEAK="${from_dir}/${cl}_peaks.narrowPeak"
if [ ! -f ${POOLED_PEAK} ]; then
    exit 1
fi

REP1_PEAK=$(reuse_pooled_peak "${from_dir}/${cl}.early_peaks.narrowPeak" ${POOLED_PEAK})
REP2_PEAK=$(reuse_pooled_peak "${from_dir}/${cl}.later_peaks.narrowPeak" ${POOLED_PEAK})
PSEUDO1_PEAK=$(reuse_pooled_peak "${from_dir}/${cl}.early.pseudo_peaks.narrowPeak" ${POOLED_PEAK})
PSEUDO2_PEAK=$(reuse_pooled_peak "${from_dir}/${cl}.later.pseudo_peaks.narrowPeak" ${POOLED_PEAK})

PooledInRep1AndRep2="${out_dir}/${cl}.PooledInRep1AndRep2.narrowPeak.gz"
PooledInPsRep1AndPsRep2="${out_dir}/${cl}.PooledInPsRep1AndPsRep2.narrowPeak.gz"
naivePeakList="${out_dir}/${cl}.naivePeakList.narrowPeak.gz"

# npeak2cl=${out_dir}/mba.whole.L4.npeak4anno.txt
# summit_list=${out_dir}/mba.whole.naiveSummitList.list 

### Naive overlap
# Find pooled peaks that overlap Rep1 and Rep2
# where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
bedtools intersect -wo -a ${POOLED_PEAK} -b ${REP1_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | bedtools intersect -wo -a stdin -b ${REP2_PEAK} \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInRep1AndRep2}

# Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2
# where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
bedtools intersect -wo -a ${POOLED_PEAK} -b ${PSEUDO1_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq \
    | bedtools intersect -wo -a stdin -b ${PSEUDO2_PEAK} -nonamecheck \
    | awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
    | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInPsRep1AndPsRep2}

# Combine peak lists
zcat ${PooledInRep1AndRep2} ${PooledInPsRep1AndPsRep2}  \
    | sort -k1,1 -k2,2n | uniq \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[0-9XY]+(?!_)' | gzip -c  > ${naivePeakList}

# Get summit 
POOLED_SUMMIT="${from_dir}/${cl}.summits.bed"
cat ${POOLED_PEAK} | awk 'BEGIN{FS=OFS="\t"}{s1=$2+$10; s2=$2+$10+1}{print $1,s1,s2,$4,$9}' > ${POOLED_SUMMIT}

naiveSummitList="${out_dir}/${cl}.naiveSummitList.bed"
join -1 1 -2 4 <(zcat ${naivePeakList} | cut -f 4 | sort) <(sort -k4,4 ${POOLED_SUMMIT}) -t$'\t' \
| awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' | sort -k1,1 -k2,2n > ${naiveSummitList} 

# if [ -f ${naiveSummitList} ] && [ -s ${naiveSummitList} ]; then
#     echo -e "${cl}\t${naiveSummitList}" >> ${summit_list}
# 	# summary peaks in cell type
# 	npeak=`cat ${naiveSummitList} | wc -l`
# 	echo -e "${cl}\t${npeak}" >> $npeak2cl
# fi

echo "${cl} done."
