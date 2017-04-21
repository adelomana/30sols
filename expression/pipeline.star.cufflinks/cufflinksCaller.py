import sys,os

def cuffdiffCaller():

    '''
    this function calls cuffdiff using specific conditions
    '''

    # working on tp1 vs tp4
    outputDir=cufflinksDir+'tp4.vs.tp1'
    controlFiles=[element for element in abundanceFiles if 'tp.1' in element and 'trna' in element]
    sampleFiles=[element for element in abundanceFiles if 'tp.4' in element and 'trna' in element]

    term1='time cuffdiff %s '%(gtfFile)
    term2=','.join(controlFiles)+' '+','.join(sampleFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    cmd=term1+term2+term3+term4+term5
    
    print('')
    print(cmd)
    print('')

    os.system(cmd)

    # working on tp1 vs tp3
    outputDir=cufflinksDir+'tp3.vs.tp1'
    controlFiles=[element for element in abundanceFiles if 'tp.1' in element and 'trna' in element]
    sampleFiles=[element for element in abundanceFiles if 'tp.3' in element and 'trna' in element]

    term1='time cuffdiff %s '%(gtfFile)
    term2=','.join(controlFiles)+' '+','.join(sampleFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    cmd=term1+term2+term3+term4+term5
    
    print('')
    print(cmd)
    print('')

    os.system(cmd)

    # working on tp1 vs tp2
    outputDir=cufflinksDir+'tp2.vs.tp1'
    controlFiles=[element for element in abundanceFiles if 'tp.1' in element and 'trna' in element]
    sampleFiles=[element for element in abundanceFiles if 'tp.2' in element and 'trna' in element]

    term1='time cuffdiff %s '%(gtfFile)
    term2=','.join(controlFiles)+' '+','.join(sampleFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    cmd=term1+term2+term3+term4+term5
    
    print('')
    print(cmd)
    print('')

    os.system(cmd)

    return None

def cuffnormCaller():

    '''
    this function calls cuffnorm using all generated abundance binary files by cuffquant
    '''

    outputDir=cufflinksDir+'allSamples'

    term1='time cuffnorm %s '%(gtfFile)
    term2=' '.join(abundanceFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    term6='-L '+','.join(labels)+' '
    cmd=term1+term2+term3+term4+term5+term6
    
    print('')
    print(cmd)
    print('')

    os.system(cmd)

    return None

def cuffquantCaller(inputFile):

    '''
    this function calls the different steps of the cufflinks pipeline.
    '''

    label=inputFile.split('/')[-2]
    outputDir=cufflinksDir+label
    
    # cuffquant
    term1='time cuffquant %s '%(gtfFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='-M %s '%maskFile
    term5='--library-type fr-firststrand '
    term6='--multi-read-correct '
    cmd=term1+term2+term3+term4+term5+term6+inputFile

    print('')
    print(cmd)
    print('')

    os.system(cmd)

    return None

# 0. defining input files
bamFilesDir='/proj/omics4tb/alomana/projects/TLR/data/BAM/'
cufflinksDir='/proj/omics4tb/alomana/projects/TLR/data/cufflinks/'
gtfFile='/proj/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'
maskFile='/proj/omics4tb/alomana/projects/TLR/data/genome/mask.gff3'
numberOfThreads=16

# 1. defining the BAM and abundance files
roots=os.listdir(bamFilesDir)
bamFiles=[bamFilesDir+element+'/Aligned.sortedByCoord.out.bam' for element in roots]
abundanceFiles=[cufflinksDir+element+'/abundances.cxb' for element in roots]
labels=[element.split('_')[-1] for element in roots]

# 2. calling cuffquantCaller 
#print('calling cuffquant...')
#for inputFile in bamFiles:
#    cuffquantCaller(inputFile)

# 3. calling cuffnorm
#print('calling cuffnorm...')
#cuffnormCaller()

# 4. calling cuffdiff
print('calling cuffdiff...')
cuffdiffCaller()
