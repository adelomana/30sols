###
### This script compares ribo-pt to RNA-seq or Ribo-seq
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    this function creates a dictionary for expression values as
    expression[trna/rbf][ribo-pt gene name][timepoint][replicate]=value
    '''

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    with open(expressionDataFile,'r') as f:

        firstLine=f.readline()
        header=firstLine.split(',')
        sampleNames=header[1:]
        sampleNames[-1]=sampleNames[-1].replace('\n','')

        for line in f:
            vector=line.split(',')

            # geneName
            geneName=vector[0].replace('_','')
            if geneName in riboPtNames:

                for i in range(len(sampleNames)):

                    # sampleType
                    sampleType=sampleNames[i].split('.')[0]
                    if sampleType not in sampleTypes:
                        sampleTypes.append(sampleType)

                    # timepoint
                    timepoint='tp.{}'.format(int(sampleNames[i].split('.')[-1]))
                    if timepoint not in timepoints:
                        timepoints.append(timepoint)

                    # replicate
                    replicate='rep.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
                    if replicate not in replicates:
                        replicates.append(replicate)

                    # value
                    value=float(vector[i+1])

                    # make sure keys exist
                    if sampleType not in expression.keys():
                        expression[sampleType]={}
                    if geneName not in expression[sampleType].keys():
                        expression[sampleType][geneName]={}
                    if timepoint not in expression[sampleType][geneName].keys():
                        expression[sampleType][geneName][timepoint]={}

                    expression[sampleType][geneName][timepoint][replicate]=value

    # sort variables
    sampleTypes.sort()
    geneNames.sort()
    timepoints.sort()
    replicates.sort()
    
    return expression,sampleTypes,timepoints,replicates

def proteinDataReader():

    '''
    This function reads data available and outputs the defined dictionary.
    '''

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC

    conditions=[]
    geneNames=[]
    timepoints=[]
    replicates=[]
    
    allFiles=os.listdir(proteinDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteinDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in conditions:
            conditions.append(condition)
        if replicate not in replicates:
            replicates.append(replicate)

        if condition not in data.keys():
            data[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}

        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')

                geneName=vector[0]
                if geneName not in geneNames:
                    geneNames.append(geneName)

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c

    # sort
    geneNames.sort()
    conditions.sort()
    replicates.sort()
    
    return data,geneNames,conditions,replicates,timepoints

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])

    riboPtNames.sort()
            
    return riboPtNames

###
### MAIN
###

# 0. user defined variables
proteinDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

# 1. read data
print('reading data...')

# 1.1. read ribo-pt names
riboPtNames=riboPtNamesReader()

# 1.2. read proteome quantification
proteinData,geneNames,fractions,replicatesProtein,timepointsProtein=proteinDataReader()

# 1.3. read transcriptome data
expression,sampleTypes,timepointsTranscript,replicatesTranscript=expressionReader()

# 2. analysis
print('performing analysis...')

timepointLate=timepointsTranscript[-1]
timepointEarly=timepointsTranscript[0]

for fraction in fractions:
    for sampleType in sampleTypes:

        print('working with {} and {}...'.format(fraction,sampleType))

        figureFileName='figure.{}.{}.pdf'.format(fraction,sampleType)

        x=[]; y=[]
        
        for name in riboPtNames:
            
            # 2.1. compute x-axis, transcript

            # compute averages for transcript late
            a=numpy.mean([expression[sampleType][name][timepointLate][replicate] for replicate in replicatesTranscript])

            # compute averages for transcript late
            b=numpy.mean([expression[sampleType][name][timepointEarly][replicate] for replicate in replicatesTranscript])

            # compute log2 fold-change
            log2fcx=a-b
            #x.append(log2fcx)

            # 2.2. compute y-axis, protein
            log2fcy=None
            try:
                log2fcy=numpy.median([proteinData[fraction][replicate]['tp4vs1'][name] for replicate in replicatesProtein])
            except:
                print('Excluding {} for no value at protein level...'.format(name))
                pass

            # 2.3 append if valid values are found
            if log2fcy != None:
                x.append(log2fcx); y.append(log2fcy)

        # 2.3. plot
        matplotlib.pyplot.plot(x,y,'o',alpha=0.5,mew=0,ms=8,color='black')
        matplotlib.pyplot.text(x,y,name)
        
        # labels
        matplotlib.pyplot.xlabel('Transcript ({}), log$_2$ FC'.format(sampleType))
        matplotlib.pyplot.ylabel('Protein ({}), log$_2$ FC'.format(fraction))

        # closing figure
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.axes().set_aspect('equal')
        matplotlib.pyplot.savefig(figureFileName)
        matplotlib.pyplot.clf()  
