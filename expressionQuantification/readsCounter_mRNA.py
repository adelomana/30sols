import sys,os

'''
this script quantifies the read counts using HTSeq. To be run in osiris.
'''

def htseqCounter(folder):

    '''
    this function calls HTSeq
    '''

    htseqExecutable='htseq-count -m union -r pos -f bam -t gene -s reverse -i locus_tag '
    bamFile=bamFilesDir+folder+'/Aligned.sortedByCoord.out.bam '
    genomeAnnotationFile='/Users/arjun/data/ncbi/halo.genomic.archive.gff'
    outputDirection=' > %s%s.txt'%(countsDir,folder)
    
    cmd=htseqExecutable+bamFile+genomeAnnotationFile+outputDirection

    print
    print cmd
    print
    os.system(cmd)

    return None

# 0. defining user variables
bamFilesDir='/Users/arjun/data/bamFiles/Total_RNA/'
countsDir='/Users/arjun/data/htseqCounts/vng_total_RNA/'

# 1. defining the BAM files
bamRoots=os.listdir(bamFilesDir)

# 2. calling HTSeq
for folder in bamRoots:
    htseqCounter(folder)
