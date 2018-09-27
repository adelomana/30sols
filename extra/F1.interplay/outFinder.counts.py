import sys
import matplotlib,matplotlib.pyplot
import numpy,numpy.linalg

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    This function creates a dictionary for expression values as
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
            if geneName not in geneNames:
                geneNames.append(geneName)

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

    return expression,sampleTypes,geneNames,timepoints,replicates

###
### MAING
###

# 0. user defined variables
theColor={}
theColor['tp.1']='red'
theColor['tp.2']='orange'
theColor['tp.3']='green'
theColor['tp.4']='blue'

# 0.1. paths
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
expression,expressionSampleTypes,expressionGeneNames,expressionTimepoints,expressionReplicates=expressionReader()

# 2. computing the figure
print('computing the analysis...')

# 2.1. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. remove outliers
sat=[]

for timepoint in expressionTimepoints:

    figureName='figures/outfinder.counts.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    
    for name in expressionGeneNames:

        # check consistency of mRNA
        mRNACounts=[]; footprintCounts=[]
        for replicate in expressionReplicates:
                      
            mRNACounts.append(expression['trna'][name][timepoint][replicate])
            footprintCounts.append(expression['rbf'][name][timepoint][replicate])
            
        # data transformations 
        m=numpy.median(mRNACounts); f=numpy.median(footprintCounts)
        if m < 0:
            m=0
        if f < 0:
            f=0
            
        r=f-m

        if f == 0:
            hollowx.append(m); hollowy.append(r)
        else:
            setx.append(m); sety.append(r)

    # performing regression analysis
    A=numpy.vstack([setx,numpy.ones(len(setx))]).T
    solution,residuals,rank,s=numpy.linalg.lstsq(A,sety)

    m=solution[0]
    c=solution[1]
    expected=list(m*numpy.array(setx)+c)

    print(timepoint,m,c)

    # computed from Matt Wall on log2
    satx=2**(numpy.array(setx))-1
    saty=(2**c)*((satx+1)**(1+m))-1

    sat.append([list(satx),list(saty)])
        
    # plotting
    matplotlib.pyplot.plot(setx,sety,'o',alpha=0.1,mew=0,color='black')
    matplotlib.pyplot.plot(hollowx,hollowy,'o',alpha=0.1,mew=0,color='tab:brown')
    
    matplotlib.pyplot.plot(setx,expected,'-',lw=2,color=theColor[timepoint])
    
    matplotlib.pyplot.xlabel('mRNA [log$_{2}$ counts]')
    matplotlib.pyplot.ylabel('RF/mRNA [log$_{2}$ ratio]')

    #matplotlib.pyplot.xlim([-0.1,20])
    #matplotlib.pyplot.ylim([-5,13])

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

# plot the model
figureName='figures/saturation.counts.pdf'
theColors=['red','orange','green','blue']

for i in range(len(sat)):
    px=sat[i][0]
    py=sat[i][1]
    matplotlib.pyplot.plot(px,py,'o',color=theColors[i])

matplotlib.pyplot.xlabel('mRNA [count]')
matplotlib.pyplot.ylabel('Predicted RF [count]')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()
