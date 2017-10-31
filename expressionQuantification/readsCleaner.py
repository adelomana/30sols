### this script calls Trimmomatic to clean the reads
import os,sys

def trimmomaticCaller(instance):
    '''
    This function deals with the trimming of the reads using Trimmomatic. Recommended options, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    '''

    print 'working with file',instance

    logFile=logFilesDir+instance+'.messagesFromTrimming.txt'

    inputFile1=rawReadsDir+instance+'_R1_001.fastq.gz'

    outputFile1=cleanReadsDir+instance+'_R1_001.clean.fastq.gz'

    cmd='java -jar /Users/arjun/Applications/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 4 -phred33 -trimlog %s %s %s ILLUMINACLIP:%s:2:30:7 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:7 MINLEN:15'%(logFile,inputFile1,outputFile1,path2Adapter)  

    print cmd
    os.system(cmd)
    print
    
    return None

# 0. defining user variables
rawReadsDir='/Users/arjun/Documents/basespace/20160421_ribosomal_frag_RNA/'
cleanReadsDir='/Users/arjun/data/cleanReads/RPF_all/'
logFilesDir='/Users/arjun/data/logFilesTrimmomatic/all_samples/'

path2Adapter='/Users/arjun/Applications/Trimmomatic-0.35/adapters/TruSeq3-SE.fa'

# 1. reading the files
foundFiles=os.listdir(rawReadsDir)
readsFiles=[element for element in foundFiles if not element.startswith('.')]

allFiles=[]
for readsFile in readsFiles:
    label=readsFile.split('/')[-1].split('_R1')[0]
    allFiles.append(label)
sequencingInstances=list(set(allFiles))

# 2. calling Trimmomatic
for instance in sequencingInstances:
    trimmomaticCaller(instance)
