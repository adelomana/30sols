#!/bin/bash 

time prefetch SRR6761663 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761664 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761665 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761666 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761667 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761668 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761669 -v --output-directory /Users/alomana/scratch/
time prefetch SRR6761670 -v --output-directory /Users/alomana/scratch/

time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761663.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761664.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761665.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761666.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761667.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761668.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761669.sra --outdir /Users/alomana/scratch/
time fastq-dump -v --split-files -W /Users/alomana/scratch/SRR6761670.sra --outdir /Users/alomana/scratch/
