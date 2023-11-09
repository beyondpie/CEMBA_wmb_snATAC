#!/bin/bash

# echo $1
# echo $2
# echo $3
# echo done

if [ ! -f $1 ]; then
    echo "$1 does not exist."
    exit 1
fi

if [ ! -f $2 ]; then
    echo "$2 does not exist."
    exit 1
fi

echo "intersect merge peak for $3 ."

intersectBed -wa -a $1 -b <(sed -e "1d" $2) -nonamecheck \
   | sort -k1,1 -k2,2n | uniq > $3

echo "done."
