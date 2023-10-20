#!/bin/bash
union_peakset=$1
if [ ! -f ${union_peakset} ];then
    echo "${union_peakset} does not exist."
    exit 1
fi

sed -e "1d " ${union_peakset} | cut -f 1-3 > $2

echo "done."
