import os,sys

# 0. user defined variables
resultsFolder='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

# 1. build a dictionary of new to old annotation
lexicon={}
with open(annotationFile,'r') as f:
    next(f)
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        info=v[-1]
        info=info.replace('\n','')
        new=None; old=None
        if 'old_locus_tag' in info:
            new=info.split('ID=')[1].split(';')[0]
            old=info.split('old_locus_tag=')[1].split('%')[0].split(';')[0]
            lexicon[new]=old

# 2. locate input files
allFiles=os.listdir(resultsFolder)
workingFiles=[element for element in allFiles if 'tp.4' in element and element[0] != '.' and 'old.annotation' not in element]

# 3. convert into old annotation
for inputFile in workingFiles:
    inputPath=resultsFolder+inputFile
    outputPath=inputPath.replace('.csv','.old.annotation.csv')

    g=open(outputPath,'w')
    f=open(inputPath,'r')

    next(f)
    for line in f:
        newName=line.split(',')[0].replace('"','')
        if newName in lexicon:
            oldName=lexicon[newName]
            newLine=line.replace(newName,oldName)
            g.write(newLine)

    f.close()
    g.close()
