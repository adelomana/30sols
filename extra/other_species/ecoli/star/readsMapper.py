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

def STARcalling(inputFile):

    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+inputFile+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    fastaFile=readsFilesDir+inputFile+'_clean.fastq'

    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s'%fastaFile   
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 1109973778'

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
readsFilesDir='/Users/alomana/scratch/clean_fastq/'
bamFilesDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/bam/'
STARexecutable='/Users/alomana/software/STAR-2.7.3a/bin/MacOSX_x86_64/STAR'
numberOfThreads=8
genomeIndexDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/genome/STARindex'
genomeFastaFile='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa'             
genomeAnnotationFile='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/genome/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3'   


# 1. recover the clean FASTQ files
print('reading FASTQ files...')
allTags=[]
allFiles=os.listdir(readsFilesDir)
for element in allFiles:
    tag=element.split('_clean.fastq')[0]
    allTags.append(tag)
inputFiles=list(set(allTags))
inputFiles.sort()

print(inputFiles)

# 2. making genome indexes
print('making genome index...')
genomeIndexer()
       
# 3. calling STAR
print('\ncalling STAR...')
for inputFile in inputFiles:
    STARcalling(inputFile)
