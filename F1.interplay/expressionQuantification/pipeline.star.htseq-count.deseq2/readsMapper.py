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
    flag6=' --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 75 --genomeSAindexNbases 8'

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

    fastaFile=readsFilesDir+inputFile+'.clean.fastq'

    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s'%fastaFile   
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5357465103'

    cmd='time '+STARexecutable+flag1+flag2+flag3+flag4+flag5
    
    print()
    print(cmd)
    print()
    os.system(cmd)
    
    return None

# 0. defining several input/output paths
readsFilesDir='/proj/omics4tb/alomana/projects/TLR/data/cleanFASTQ/'
bamFilesDir='/proj/omics4tb/alomana/projects/TLR/data/BAM/'
STARexecutable='/proj/omics4tb/alomana/software/STAR-2.5.4b/bin/Linux_x86_64/STAR'
genomeIndexDir='/proj/omics4tb/alomana/projects/TLR/data/genomeIndex'
genomeFastaFile='/proj/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.fasta'              
genomeAnnotationFile='/proj/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'   
numberOfThreads=16

# 1. recover the clean FASTQ files
print('reading FASTQ files...')
allTags=[]
allFiles=os.listdir(readsFilesDir)
for element in allFiles:
    tag=element.split('.clean.fastq')[0]
    allTags.append(tag)
inputFiles=list(set(allTags))
inputFiles.sort()

# 2. making genome indexes
print('making genome index...')
genomeIndexer()
       
# 3. calling STAR
print('calling STAR...')
for inputFile in inputFiles:
    STARcalling(inputFile)
