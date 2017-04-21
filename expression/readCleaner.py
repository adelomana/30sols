###
### this script calls Trimmomatic to clean raw reads
###

import os,sys,time

def trimmomaticCaller(instance):
    
    '''
    this function deals with the trimming of the reads using Trimmomatic. Recommended options, ILLUMINACLIP:path2AdaptersFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    '''

    print('working with file {}'.format(instance))

    logFile=logFilesDir+instance+'.messagesFromTrimming.txt'
    inputFile=rawReadsDir+instance+'.fastq'
    outputFile=cleanReadsDir+instance+'.clean.fastq'

    if 'total' in instance:
        minLen=36
    elif 'ribosomal' in instance:
        minLen=15
    else:
        raise RuntimeError("No flag found.")

    cmd='time {} -jar {} SE -threads 4 -phred33 -trimlog {} {} {} ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{}'.format(javaPath,trimmomaticPath,logFile,inputFile,outputFile,path2Adapter,minLen)

    print(cmd)
    os.system(cmd)
    print()
            
    return None

### MAIN

# 0. defining user variables
rawReadsDir='/Volumes/omics4tb/alomana/projects/TLR/data/FASTQ/'
cleanReadsDir='/Volumes/omics4tb/alomana/projects/TLR/data/cleanFASTQ3/'
logFilesDir='/Users/alomana/scratch/trimmomaticLogs/'
path2Adapter='/Users/alomana/software/Trimmomatic-0.36/adapters/TruSeq3-SE.fa'

javaPath='/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java'
trimmomaticPath='/Users/alomana/software/Trimmomatic-0.36/trimmomatic-0.36.jar'

# 1. reading the files
tag='.fastq'
fastqFiles=os.listdir(rawReadsDir)
readFiles=[element for element in fastqFiles if tag in element]

# 2. calling Trimmomatic
for readFile in readFiles:
    instance=readFile.split(tag)[0]
    trimmomaticCaller(instance)
