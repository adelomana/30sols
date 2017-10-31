import os,sys

def STARcalling(inputFile):
    
    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+inputFile+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    readF1=readsFilesDir+inputFile+'_L001_R1_001.clean.fastq.gz'
    readF2=readsFilesDir+inputFile+'_L002_R1_001.clean.fastq.gz'
    readF3=readsFilesDir+inputFile+'_L003_R1_001.clean.fastq.gz'
    readF4=readsFilesDir+inputFile+'_L004_R1_001.clean.fastq.gz'
    
    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesCommand gunzip -c --readFilesIn %s,%s,%s,%s'%(readF1,readF2,readF3,readF4) 
    flag4=' --outFileNamePrefix %s'%finalDir
    flag7=' --seedSearchStartLmax 13 --winAnchorMultimapNmax 200 --outFilterMatchNmin 15 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate'
    
    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag7
    
    print
    print cmd
    print
    os.system(cmd)

    return None

# 0. defining several input/output paths
readsFilesDir='/Users/arjun/data/cleanReads/RPF_all/'
bamFilesDir='/Users/arjun/data/bamFiles/RPF_all_min15/'
STARexecutable='/Users/arjun/STAR/bin/MacOSX_x86_64/STAR'
genomeIndexDir='/Users/arjun/data/genomeIndex'
numberOfThreads=4

# 1. recover the clean FASTQ files
print 'reading FASTQ files...'
allTags=[]
allFiles=os.listdir(readsFilesDir)
for element in allFiles:
    tag=element.split('_L00')[0]
    allTags.append(tag)
inputFiles=list(set(allTags))

# 2. calling STAR
print 'calling STAR...'
for inputFile in inputFiles:
    STARcalling(inputFile)
