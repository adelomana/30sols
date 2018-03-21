#!/bin/bash

#$ -N profiler
#$ -o /proj/omics4tb/alomana/scratch/messages.o.txt
#$ -e /proj/omics4tb/alomana/scratch/messages.e.txt
#$ -P Bal_alomana
#$ -pe serial 58
#$ -q baliga
#$ -S /bin/bash

cd /users/alomana
source .bash_profile

time python /proj/omics4tb/alomana/projects/TLR/src/deployment/profiler.py
