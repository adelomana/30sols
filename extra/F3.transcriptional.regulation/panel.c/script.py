###
### this script creates a figure of FC and average expression for RNA-seq and Ribo-seq data sets
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def annotationReader():

    '''
    This function reads the annotation of microbes online and return a diccionary with GI identifiers to be used in DAVID 6.7 (these IDs are not recognized by DAVID 6.8).
    '''

    annotationMap={}

    with open(annotationFile,'r') as f:
        next(f)
        for line in f:
            v=line.split('\t')
            GI=v[2]
            sysName=v[7]
            annotationMap[sysName]=GI

    return annotationMap

def expressionReader():

    '''
    This function creates a dictionary for expression values as:
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
            geneName=vector[0]
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

def grapher():

    '''
    This function builds the figure.
    '''

    timepointLate=timepoints[-1]
    timepointEarly=timepoints[0]
        
    rnax=[]; rnay=[]
    ribox=[]; riboy=[]

    hdg_FC=[]; hdg_E=[]    # highly-downregulated group
    dg_FC=[]; dg_E=[]      # downregulated group

    for geneName in riboPtNames:

        # f.1. computing the FCs

        # compute averages for RNA-seq late
        x=numpy.mean([expression['trna'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for RNA-seq early
        y=numpy.mean([expression['trna'][geneName][timepointEarly][replicate] for replicate in replicates])

        # compute averages for Ribo-seq late
        z=numpy.mean([expression['rbf'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for Ribo-seq early
        w=numpy.mean([expression['rbf'][geneName][timepointEarly][replicate] for replicate in replicates])

        # compute fold-changes
        rnaFC=x-y
        riboFC=z-w

        # append variables for plotting
        rnay.append(rnaFC)
        riboy.append(riboFC)

        # f.2. computing the average values
        x=[]; y=[]
        for timepoint in timepoints:
            x.append(numpy.mean([expression['trna'][geneName][timepoint][replicate] for replicate in replicates]))
            y.append(numpy.mean([expression['rbf'][geneName][timepoint][replicate] for replicate in replicates]))

        # making a log10 transformation of log2 normalized counts
        rnax.append(numpy.log10(2**numpy.median(x)))
        ribox.append(numpy.log10(2**numpy.median(y)))

        # f.3. descriptive statistics for two groups
        if rnaFC < -3:
            hdg_FC.append(rnaFC)
            hdg_E.append(numpy.log10(2**numpy.median(x)))
            print('group B','\t',riboPtMap[geneName],'\t',annotationMap[riboPtMap[geneName]],'\t',rnaFC)
        else:
            dg_FC.append(rnaFC)
            dg_E.append(numpy.log10(2**numpy.median(x)))
            print('group A','\t',riboPtMap[geneName],'\t',annotationMap[riboPtMap[geneName]],'\t',rnaFC)

    # printing the statistics            
    print('average HDG:',len(hdg_FC),numpy.median(hdg_FC),numpy.median(hdg_E),10**numpy.median(hdg_E))
    print('average DG:',len(dg_FC),numpy.median(dg_FC),numpy.median(dg_E),10**numpy.median(dg_E))

    # making the figure
        
    theSize=8

    f, (ax1, ax2) = matplotlib.pyplot.subplots(1, 2, sharey=True)
    ax1.plot(rnax,rnay,'o',alpha=0.5,mew=0,ms=theSize,color='black',label='RNA-seq')    
    ax1.set_xlim([1.5,4.5])
    ax1.set_ylim([-5,0])
    ax1.set_xlabel('Median over time,')
    ax1.set_ylabel('log$_2$ FC')
    ax1.set_xticks([2,3,4])
    ax1.grid(alpha=0.5, ls=':')
    ax1.legend(framealpha=1,loc=1,ncol=1,fontsize=14)
    
    ax2.plot(ribox,riboy,'o',alpha=0.5,mew=0,ms=theSize,color='red',label='Ribo-seq')
    ax2.set_xlim([1.5,4.5])
    ax2.set_ylim([-5,0])
    ax2.set_xlabel('log$_2$(normalized counts)')
    ax2.set_xticks([2,3,4])
    ax2.grid(alpha=0.5, ls=':')
    ax2.legend(framealpha=1,loc=1,ncol=1,fontsize=14)
    
    figureName='figure.changes.pdf'
    f.tight_layout()
    f.savefig(figureName)
    f.clf()
    
    return None

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    riboPtMap={}
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])
            riboPtMap[vector[0]]=vector[1]
            
    return riboPtNames,riboPtMap

### MAIN

# 0. user defined variables
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/genome.annotation.txt'

# 1. read the data
print('reading data...')
riboPtNames,riboPtMap=riboPtNamesReader()
expression,sampleTypes,timepoints,replicates=expressionReader()

annotationMap=annotationReader()

# 2. plot the data
print('plotting data...')
grapher()

print('... done.')
