#!/bin/bash

#Author: Yang Li <yal054@health.ucsd.edu>
#File: run.sh
#Create Date: Sat Feb 26 10:26:52 PST 2022

snakemake -p --rerun-incomplete -k -j 128 --cluster "qsub -l {cluster.ppn} -l {cluster.time} -N {params.jobname} -q {cluster.queue} -o pbslog/{params.jobname}.pbs.out -e pbslog/{params.jobname}.pbs.err" --jobscript jobscript.pbs --jobname "{rulename}.{jobid}.pbs" --cluster-config cluster.json 2>run.log
