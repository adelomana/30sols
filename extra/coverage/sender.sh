#!/bin/bash

#$ -N profiler
#$ -o /proj/omics4tb/alomana/scratch/messages.o.txt
#$ -e /proj/omics4tb/alomana/scratch/messages.e.txt
#$ -pe smp 25
#$ -S /bin/bash

cd /users/alomana
source .bash_profile

time python /proj/omics4tb/alomana/projects/TLR/src/deployment/profiler.py
