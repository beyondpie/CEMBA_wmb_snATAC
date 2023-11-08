#!/bin/bash

usage() {
    echo "Usage: $0 -c [cicero_dir] -m [meta_col]"
    exit 1
}

while getopts ":c:a:g:m:" o; do
    case "${o}" in
        c)
            cicero_dir=${OPTARG};;
        m)
            meta_col=${OPTARG};;
        *)
            usage ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${cicero_dir}" ] || [ -z "${meta_col}" ]; then
    usage
fi

echo "Merge distal-proximal conns: ${meta_col}"
out_all_conns_sta="${cicero_dir}/mba.whole.${meta_col}.merge.conns.sta.all"
cat ${cicero_dir}/${meta_col}.*.conns.sta \
    | sed -e "1i cluster\ttotConns\tddConns\tppConns\tpdConns\tgeneN\tcreN" \
          > ${out_all_conns_sta}
out_all_pdc_sta="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.sta.all"
cat ${cicero_dir}/${meta_col}.*.pdc.sta \
    | sed -e "1i cluster\ttotc\tgeneN\tcreN" \
          > ${out_all_pdc_sta}

## All pdc without any filtering
out_all_pdc="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.all"
cat ${cicero_dir}/${meta_col}.*.fitConns.res.alignv1.pdc \
    | grep -v 'anno1' \
    | sed -e "1i peak1\tcre1\tanno1\tpeak2\tcre2\tanno2\tcoaccess\tpval\tfdr\tcluster\tgene" \
          > ${out_all_pdc}

# NOTE: since differen subclass may have different scores for the same
# pdc, then there are repeat pdcs in this file.
out_all_bedpe="${cicero_dir}/mba.whole.${meta_col}.merge.bedpe.all"
sed '1d' ${out_all_pdc} \
    | awk 'BEGIN{FS=OFS="\t"}{split($1,a,"_"); split($4,b,"_"); print a[1],a[2],a[3],b[1],b[2],b[3],$11"|"$5,$7,".","."}'\
    | sort -k1,1 -k2,2n | uniq > ${out_all_bedpe}

out_all_pdcpair="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.pair.all"
awk 'BEGIN{FS=OFS="\t"}{print $11"|"$5,$2}' ${out_all_pdc} \
    | sed '1d' | sort | uniq > ${out_all_pdcpair}

echo "# of cCRE"
cut -f 5 ${out_all_pdc} | sed '1d' | sort | uniq | wc -l
echo "# of genes"
cut -f 11 ${out_all_pdc} | sed '1d' | sort | uniq | wc -l
echo "# of pairs"
sed '1d' ${out_all_pdc} | sort | uniq | wc -l

echo "Merge all proximal-distal pairs' distances"
out_all_dist="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.dist.all"
cat ${cicero_dir}/${meta_col}.*.pdc.dist \
    | sed -e "1i conns\tdistance\tcluster" > ${out_all_dist}

echo "Merge all genes' stat per cluster"
out_all_peak2gene="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.peak2gene.all"
cat ${cicero_dir}/${meta_col}.*.pdc.peak2gene \
    | sed -e "1i gene\tcnt\tcluster" > ${out_all_peak2gene}

echo "Merge all peaks' stat per cluster"
out_all_gene2peak="${cicero_dir}/mba.whole.${meta_col}.merge.pdc.gene2peak.all"
cat ${cicero_dir}/${meta_col}.*.pdc.gene2peak \
    | sed -e "1i cCREs\tcnt\tcluster" > ${out_all_gene2peak}

