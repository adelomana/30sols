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
    flag3='-t mRNA'
    flag5='-i ID'    
        
    #flag4='-s yes'
    flag4='-s no'
    
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
bamFilesDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/bam/'
countsDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/counts/'
genomeAnnotationFile='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3' 
numberOfThreads=6

# 1. defining the BAM files
samples=os.listdir(bamFilesDir)

print(samples)

# 2. calling HTSeq in a parallel environment
#for sample in samples:
#    htseqCounter(sample)
#    sys.exit()
    
hydra=multiprocessing.pool.Pool(numberOfThreads)
hydra.map(htseqCounter,samples)
