'''
this script quantifies the read counts using HTSeq. To be run in osiris.
'''

import sys,os
import multiprocessing,multiprocessing.pool

def htseqCounter(sample):

    '''
    this function calls HTSeq
    '''

    htseqExecutable='htseq-count'
    flag1='-m union'
    flag2='-f bam'
    flag3='-t gene'
    flag5='-i ID'    
    
    if 'trna' in sample:
        flag4='-s reverse'
    elif 'rbf' in sample:
        flag4='-s yes'
    else:
        print('error a')
        sys.exit()
    
    inputFile=bamFilesDir+sample+'/Aligned.sortedByCoord.out.bam'
    outputDirection='> {}{}.txt'.format(countsDir,sample)

    cmd=' '.join(['time ',htseqExecutable,flag1,flag2,flag3,flag4,flag5,inputFile,genomeAnnotationFile,outputDirection])

    print()
    print(cmd)
    print()
    
    os.system(cmd)

    return None

###
### README
###

# 0. defining user variables
bamFilesDir='/Volumes/omics4tb/alomana/projects/TLR/data/BAM/'
countsDir='/Volumes/omics4tb/alomana/projects/TLR/data/countsArjun/'
genomeAnnotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'
numberOfThreads=4

# 1. defining the BAM files
samples=os.listdir(bamFilesDir)

# 2. calling HTSeq in a parallel environment
hydra=multiprocessing.pool.Pool(numberOfThreads)
hydra.map(htseqCounter,samples)
