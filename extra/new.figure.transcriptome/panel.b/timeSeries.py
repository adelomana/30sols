###
### this script performs analysis as Fig 1b of Schafer et al. (PMID, 26007203)
###

import sys,numpy
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

def riboPtNamesReader():

    '''
    this function reads the ribosomal protein names
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
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
theColors=['orange','green','blue']

# 1. read data
riboPtNames=riboPtNamesReader()
expression,sampleTypes,timepoints,replicates=expressionReader()

# 2. process data
values=[]
for i in range(len(timepoints)-1):

    indexLate=i+1
    timepointLate=timepoints[indexLate]
    timepointEarly=timepoints[0]
    comparison='TP {}/TP {}'.format(timepointLate[-1],timepointEarly[-1])
    print(comparison)

    # work per gene
    allFCx=[]; allFCy=[]
    for geneName in riboPtNames:

        # compute averages for RNA-seq late
        x=numpy.mean([expression['trna'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for RNA-seq early
        y=numpy.mean([expression['trna'][geneName][timepointEarly][replicate] for replicate in replicates])
        
        # compute averages for Ribo-seq late
        z=numpy.mean([expression['rbf'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for Ribo-seq early
        w=numpy.mean([expression['rbf'][geneName][timepointEarly][replicate] for replicate in replicates])
        
        # compute fold-changes
        fcx=x-y
        fcy=z-w

        if fcx != 0 or fcy != 0:
            allFCx.append(fcx); allFCy.append(fcy)
        else:
            print('excluded:\t',geneName,fcx,fcy)

    # plot values
    matplotlib.pyplot.plot(allFCx,allFCy,'o',alpha=0.5,mew=0,ms=8,color=theColors[i],label=comparison)

# closing the figure
matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

matplotlib.pyplot.legend(markerscale=1.5,framealpha=1)

matplotlib.pyplot.xlim([-5,0.5])
matplotlib.pyplot.ylim([-5,0.5])

matplotlib.pyplot.plot([-8,4],[0,0],':',color='black',alpha=0.5)
matplotlib.pyplot.plot([0,0],[-10,2],':',color='black',alpha=0.5)
matplotlib.pyplot.plot([-8,2],[-8,2],':',color='black',alpha=0.5)

figureName='figure.all.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()        
