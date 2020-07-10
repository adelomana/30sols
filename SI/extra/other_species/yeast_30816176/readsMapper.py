import os,sys

'''
This script finds the clean FASTQ files and calls STAR for the reads alignment.
'''

def genomeIndexer():

    '''
    this function creates the genome index. Should be run only once.
    '''

    flag1=' --runMode genomeGenerate'
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --genomeDir %s'%genomeIndexDir
    flag4=' --genomeFastaFiles %s'%genomeFastaFile
    flag5=' --sjdbGTFfile %s'%genomeAnnotationFile
    flag6=' --sjdbOverhang 49 --genomeSAindexNbases 8'

    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag5+flag6

    print()
    print(cmd)
    print()
    os.system(cmd)

    return None

def STARcalling(sample):

    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+sample+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    fastq_file=readsFilesDir+sample+'.fastq'

    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s'%fastq_file
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 2719138304'

    #! consider ulimit -n 512
    
    cmd='time '+STARexecutable+flag1+flag2+flag3+flag4+flag5
    
    print()
    print(cmd)
    print()
    os.system(cmd)
    
    return None

###
### MAIN
###

# 0. defining several input/output paths
readsFilesDir='/Volumes/omics4tb2/alomana/projects/TLR/data/sand/RiboSeqPy-master/4-Subtracted/'
bamFilesDir='/Volumes/omics4tb2/alomana/projects/TLR/results/yeast_358309644/bam/'
STARexecutable='/Users/alomana/software/STAR-2.7.3a/bin/MacOSX_x86_64/STAR'
numberOfThreads=8

genomeIndexDir='/Volumes/omics4tb2/alomana/projects/TLR/data/sand/annotation/starIndex'
genomeFastaFile='/Volumes/omics4tb2/alomana/projects/TLR/data/sand/annotation/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'              
genomeAnnotationFile='/Volumes/omics4tb2/alomana/projects/TLR/data/sand/annotation/Saccharomyces_cerevisiae.R64-1-1.46.gff3'   

# 1. recover the clean FASTQ files
print('reading FASTQ files...')
samples=['WTS1','WTS3','64','65']

# 2. making genome indexes
#print('making genome index...')
#genomeIndexer()
       
# 3. calling STAR
print('\ncalling STAR...')
for sample in samples:
    STARcalling(sample)
