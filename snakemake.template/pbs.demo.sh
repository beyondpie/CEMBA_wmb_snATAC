#!/bin/bash

#PBS -q [QueueName]
#PBS -N [JobName]
#PBS -l nodes=1:ppn=2,walltime=24:00:00
#PBS -V
#PBS -M [Your-email]
#PBS -m a
#PBS -A ren-group
#PBS -j oe

[Your bash script content]
