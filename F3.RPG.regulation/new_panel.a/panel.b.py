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
theColors=['orange','green','blue']
theTimePoints=[2,3,4]
sampleTypes=['trna','rbf']

# 1. read data
riboPtNames=riboPtNamesReader()
expression=expressionReader()

# 2. process data
values=[]
for i in range(len(theTimePoints)):

    comparison='tp.{}_vs_tp.1'.format(theTimePoints[i])
    comparisonLabel='TP {}/TP 1'.format(theTimePoints[i])

    # work per gene
    allFCx=[]; allFCy=[]
    for geneName in riboPtNames:

        # compute fold-changes
        fcx=expression[comparison]['trna'][geneName]
        fcy=expression[comparison]['rbf'][geneName]

        if fcx != 0 or fcy != 0:
            allFCx.append(fcx); allFCy.append(fcy)
        else:
            print('excluded:\t',geneName,fcx,fcy)

    # plot values
    matplotlib.pyplot.plot(allFCx,allFCy,'o',alpha=0.5,mew=0,ms=8,color=theColors[i],label=comparisonLabel)

# closing the figure
matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

#matplotlib.pyplot.legend(markerscale=1.5,framealpha=1)

matplotlib.pyplot.xticks([-6,-5,-4,-3,-2,-1,0])

matplotlib.pyplot.xlim([-6.5,0.5])
matplotlib.pyplot.ylim([-6.5,0.5])

#matplotlib.pyplot.plot([-6,0],[0,0],':',color='black',alpha=0.5)
#matplotlib.pyplot.plot([0,0],[-6,0],':',color='black',alpha=0.5)
matplotlib.pyplot.plot([-6,0],[-6,0],':',color='black',alpha=0.5)

figureName='figure.all.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()        
