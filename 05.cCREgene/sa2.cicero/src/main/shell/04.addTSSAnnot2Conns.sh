#!/bin/bash

usage() {
    echo "Usage: $0 -c [cicero_dir] -a [peakannot_file] -g [group] -m [meta_col]"
    exit 1
}

while getopts ":c:a:g:m:" o; do
    case "${o}" in
        c)
            cicero_dir=${OPTARG};;
        a)
            peakannot_file=${OPTARG};;
        g)
            group=${OPTARG};;
        m)
            meta_col=${OPTARG};;
        *)
            usage ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${cicero_dir}" ] || [ -z "${peakannot_file}" ]; then
    usage
fi

if [ -z "${group}" ] || [ -z "${meta_col}" ]; then
    usage
fi

echo "Distal-proximal filtering: ${group}"

fitConn_file="${cicero_dir}/cicero.${meta_col}.${group}.fitConns.res.txt"
if [ ! -f "${fitConn_file}" ]; then
    echo -e "${fitConn_file} does not exist."
    usage
fi

fitConn_alignv1_file="${cicero_dir}/${meta_col}.${group}.fitConns.res.alignv1.txt"
cut -f 1,2,3,5,6 ${fitConn_file} > ${fitConn_alignv1_file}

out_anno="${cicero_dir}/${meta_col}.${group}.fitConns.res.anno"
out_sel="${cicero_dir}/${meta_col}.${group}.fitConns.res.sel"
out_sta="${cicero_dir}/${meta_col}.${group}.conns.sta"
out_pdc="${cicero_dir}/${meta_col}.${group}.fitConns.res.pdc"
out_alignv1_pdc="${cicero_dir}/${meta_col}.${group}.fitConns.res.alignv1.pdc"
out_bedpe="${cicero_dir}/${meta_col}.${group}.pdc.bedpe"
out_pdcsta="${cicero_dir}/${meta_col}.${group}.pdc.sta"
out_peak2gene="${cicero_dir}/${meta_col}.${group}.pdc.peak2gene"
out_gene2peak="${cicero_dir}/${meta_col}.${group}.pdc.gene2peak"
out_pdcdist="${cicero_dir}/${meta_col}.${group}.pdc.dist"

echo "use tssAnnot file to annot conns"
join -1 2 -2 1 <(sort -k2,2 ${fitConn_alignv1_file}) <(sort -k1,1 ${peakannot_file}) -t$'\t' \
    | sort -k2,2 | join -1 2 -2 1 - <(sort -k1,1 ${peakannot_file}) -t$'\t' \
    | awk 'BEGIN{FS=OFS="\t"}{print $1,$9,$10,$11,$2,$6,$7,$8,$3,$4,$5}' \
    | sed -e "1i peak1\tcre1\tclass1\tgene1\tpeak2\tcre2\tclass2\tgene2\tcoaccess\tpval\tfdr" \
          > ${out_anno}

if [ ! `wc -l ${out_anno} | awk '{print $1}'` -ge "2" ]; then
    echo -e "${out_anno} is empty."
    usage
fi

# Filtering cicero based FDR.
echo "filter conns based on FDR"
awk '$11<=0.001' ${out_anno} > ${out_sel}

echo "stat distal - proximal"
tot=`cut -f 2,3,6,7 ${out_sel} | sed '1d' | sort | uniq | wc -l`
ddc=`cut -f 2,3,6,7 ${out_sel} | sed '1d' | sort | uniq | awk '($2=="distal" && $4=="distal")' | wc -l`
ppc=`cut -f 2,3,6,7 ${out_sel} | sed '1d' | sort | uniq | awk '($2=="proximal" && $4=="proximal")' | wc -l`
pdc=`cut -f 2,3,6,7 ${out_sel} | sed '1d' | sort | uniq | awk '($2!=$4)' | wc -l`

echo "distal - proximal: num of gene / num of peak"
geneN=`cut -f 2,3,4,6,7,8  ${out_sel}  | sed '1d' \
| awk '$2!=$5' | cut -f 3,6 | tr '\t' '\n' |grep -v "nan" | sort | uniq | wc -l`
creN=`cut -f 2,3,4,6,7,8  ${out_sel}  | sed '1d' | awk '$2!=$5' \
| awk 'BEGIN{FS=OFS="\t"}{print $1"|"$2,$4"|"$5}' \
| tr '\t' '\n' | grep -v "proximal" | sort | uniq | wc -l`

echo -e "${group}\t${tot}\t${ddc}\t${ppc}\t${pdc}\t${geneN}\t${creN}" > ${out_sta}


## RAW pdc without any filtering after cicero.
echo "filter distal-proximal conns."
awk 'BEGIN{FS=OFS="\t"}($3!=$7)' ${out_sel} \
    | awk 'BEGIN{FS=OFS="\t"}(NR>1){if($3=="proximal"){print $0}else{print $5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}}'\
          > ${out_pdc}

echo "align v1 pdc format"
awk -v group=${group} 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$5,$6,$7,$9,$10,$11,group,$4}' ${out_pdc} \
    | sed -e "1i peak1\tcre1\tanno1\tpeak2\tcre2\tanno2\tcoaccess\tpval\tfdr\tgroup\tgene" \
          > ${out_alignv1_pdc}

echo "get bedpe"
sed '1d' ${out_alignv1_pdc} \
    | awk 'BEGIN{FS=OFS="\t"}{split($1,a,"_"); split($4,b,"_"); print a[1],a[2],a[3],b[1],b[2],b[3],$11"|"$5,$7,".","."}' \
    | sort -k1,1 -k2,2n | uniq > ${out_bedpe}

## FIXME: out_alignv1_pdc has some repeats, should use out_bedpe
echo "summary cicero conns"
totc=`sed '1d' ${out_alignv1_pdc} | wc -l`
geneN=`cut -f 11 ${out_alignv1_pdc} | sed '1d' | sort | uniq | wc -l `
creN=`cut -f 5 ${out_alignv1_pdc} | sed '1d' | sort | uniq | wc -l`
echo -e "${group}\t${totc}\t${geneN}\t${creN}" > ${out_pdcsta}

echo "summary #peaks to gene"
sed '1d' ${out_alignv1_pdc} \
    | cut -f 11 | sort | uniq -c \
    | awk 'BEGIN{OFS="\t"}{print $2,$1}' \
    | sed -e "s/$/\t${group}/g" > ${out_peak2gene}

echo "summary #gene have peak"
sed '1d' ${out_alignv1_pdc} | cut -f 5 | sort | uniq -c \
    | awk 'BEGIN{OFS="\t"}{print $2,$1}' \
    | sed -e "s/$/\t${group}/g" > ${out_gene2peak}


echo "summary dist between peak-gene"
cat ${out_bedpe} | awk 'BEGIN{FS=OFS="\t"}
  function abs(v) {return v < 0 ? -v : v}
  {print $7, abs(($5 + ($6-$5)/2)-($2 + ($3-$2)/2))}' \
    | sed -e "s/$/\t${group}/g" > ${out_pdcdist}

echo "Done."
