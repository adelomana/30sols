###
### this script performs analysis as Fig 1b of Schafer et al. (PMID, 26007203)
###

import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    This function creates a dictionary for expression values as:
    expression[comparison][trna/rbf][ribo-pt gene name]=foldchange
    '''

    expression={}

    for timePoint in theTimePoints:
        comparison='tp.{}_vs_tp.1'.format(timePoint)
        expression[comparison]={}
        
        for sampleType in sampleTypes:
            expressionDataFile=expressionDataDir+'significance.{}.condition_{}.csv'.format(sampleType,comparison)
            expression[comparison][sampleType]={}

            with open(expressionDataFile,'r') as f:
                next(f)
                for line in f:
                    vector=line.split(',')

                    geneName=vector[0].replace('"','')
                    log2FC=float(vector[2])

                    if geneName in riboPtNames:
                        expression[comparison][sampleType][geneName]=log2FC

    return expression

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
            
    return riboPtNames

### M A I N

# 0. user defined variables
expressionDataDir='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
theColors=['blue','green','orange']
theTimePoints=[4,3,2]
sampleTypes=['trna','rbf']

# 1. read data
riboPtNames=riboPtNamesReader()
expression=expressionReader()

# 2. process data

# 2.1. process last time point to discriminate group A vs group B
index=0
comparison='tp.{}_vs_tp.1'.format(theTimePoints[index])
comparisonLabel='TP {}/TP 1'.format(theTimePoints[index])
allFCxA=[]; allFCyA=[]
allFCxB=[]; allFCyB=[]
groupA=[]
for geneName in riboPtNames:
    fcx=expression[comparison]['trna'][geneName]
    fcy=expression[comparison]['rbf'][geneName]

    if fcx != 0 or fcy != 0:
        if fcx < -4:
            groupA.append(geneName)
            allFCxA.append(fcx); allFCyA.append(fcy)
        else:
            allFCxB.append(fcx); allFCyB.append(fcy)
    else:
        print('excluded:\t',geneName,fcx,fcy)

# plot values
matplotlib.pyplot.plot(allFCxA,allFCyA,'s',alpha=0.5,mew=0,ms=8,color=theColors[index],label=comparisonLabel)
matplotlib.pyplot.plot(allFCxB,allFCyB,'D',alpha=0.5,mew=0,ms=8,color=theColors[index],label=comparisonLabel)

# 2.2. process middle time points knowing group labels
for index in [1,2]:
    comparison='tp.{}_vs_tp.1'.format(theTimePoints[index])
    comparisonLabel='TP {}/TP 1'.format(theTimePoints[index])
    allFCxA=[]; allFCyA=[]
    allFCxB=[]; allFCyB=[]

    for geneName in riboPtNames:
        fcx=expression[comparison]['trna'][geneName]
        fcy=expression[comparison]['rbf'][geneName]

        if fcx != 0 or fcy != 0:
            if geneName in groupA:
                allFCxA.append(fcx); allFCyA.append(fcy)
            else:
                allFCxB.append(fcx); allFCyB.append(fcy)
        else:
            print('excluded:\t',geneName,fcx,fcy)

    # plot values
    matplotlib.pyplot.plot(allFCxA,allFCyA,'s',alpha=0.5,mew=0,ms=8,color=theColors[index],label=comparisonLabel)
    matplotlib.pyplot.plot(allFCxB,allFCyB,'D',alpha=0.5,mew=0,ms=8,color=theColors[index],label=comparisonLabel)

# close figure
figureName='figure.all.pdf'
matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')
matplotlib.pyplot.xticks([-6,-5,-4,-3,-2,-1,0])
matplotlib.pyplot.xlim([-6.5,0.5])
matplotlib.pyplot.ylim([-6.5,0.5])
matplotlib.pyplot.plot([-6,0],[-6,0],':',color='black',alpha=0.5)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()        
