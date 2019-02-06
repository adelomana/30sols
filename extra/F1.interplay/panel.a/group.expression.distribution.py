###
### This script makes a figure of the distribution of expression for the different regulatory groups
###

import os,sys,numpy
import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':24,'font.family':'Arial','xtick.labelsize':18,'ytick.labelsize':18})

def histogrammer(theData):

    '''
    This function creates a histogram.
    '''    

    x=[]; y=[]
    
    binSize=0.1
    left=0
    right=5
    rightBins=numpy.arange(left+binSize,right+binSize,binSize)
    n,bins=numpy.histogram(theData,bins=rightBins)

    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        x.append(center)
    x.pop()
    y=numpy.array(n)
    y=list(y/float(sum(y)))

    return x,y

def transcriptomicsReader():

    '''
    This function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    geneNames=[]; timepoints=[]; replicates=[]
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

            if replicate not in replicates:
                replicates.append(replicate)
            if timepoint not in timepoints:
                timepoints.append(timepoint)

            if fraction not in data.keys():
                data[fraction]={}
            if replicate not in data[fraction].keys():
                data[fraction][replicate]={}
            if timepoint not in data[fraction][replicate].keys():
                data[fraction][replicate][timepoint]={}
            
        for line in f:
            vector=line.split('\t')[:-1]
            values=[float(element) for element in vector[1:]]
            geneName=vector[0].replace('_','')
            if geneName not in geneNames:
                geneNames.append(geneName)
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data,geneNames,timepoints,replicates

###
### MAIN
###

# 0. user defined variables

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA expression data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 1.2. reading group membership
geneSets={}
elements=os.listdir('results')
for element in elements:
    tag=element.split('results.')[1].split('.txt')[0]
    geneSets[tag]=[]

    # read file
    file2read='results/{}'.format(element)
    with open(file2read,'r') as f:
        for line in f:
            v=line.split('\t')
            if v[2] not in ['gene-VNGRS13150','gene-VNGRS00005','gene-VNGRS09790','gene-VNGRS10040','gene-VNGRS09800','gene-VNGRS09805','gene-VNGRS03925']:
                geneSets[tag].append(v[2])

# add manuall a group with all genes
longNames=['gene-'+element for element in geneNames]
geneSets['all']=longNames
            
# 2. convert group memmberships into expression distributions
print('converting group memberships into expression distributions...')
expressionDistributions={}
for element in geneSets.keys():
    expressionDistributions[element]=[]

    # check consistency of mRNA
    count=0
    for geneName in geneSets[element]:
        mRNA_TPMs=[]
        shortGeneName=geneName.split('gene-')[1]
        for replicate in replicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate]['tp.1'][shortGeneName])
            
        # data transformations and quality check
        log10M=numpy.log10(numpy.array(mRNA_TPMs)+1)
        log2M=numpy.log2(numpy.array(mRNA_TPMs)+1)

        # noise
        if numpy.max(log2M) > numpy.log2(10+1): # if expression is below 1 TPMs, don't consider noise
            sem=numpy.std(log2M)/numpy.sqrt(len(log2M))
            rsem_mRNA=sem/numpy.mean(log2M)
        else:
            rsem_mRNA=0

        if rsem_mRNA < 0.3:
            m=numpy.median(log10M)
            expressionDistributions[element].append(m)
            
# 3. define significance of deviation
print('running hypothesis test of deviation...')

groupLabels=list(expressionDistributions.keys())
groupLabels.sort()
groupLabels.remove('dubious')

theColors=['black','black','black','blue','green','orange','red','yellow']
theLineStyle=['-',':','--','-','-','-','-','-']

# run specific groups

#groupLabels=['orange','green']
#theColors=['orange','green']

#groupLabels=['black.minus', 'black.plus','blue','red']
#theColors=['black','black','blue','red']

groupLabels=['yellow']
theColors=['yellow']

# make a figure of the overal distribution
x,y=histogrammer(expressionDistributions['all'])
matplotlib.pyplot.plot(x,y,'-',color='black',lw=1)
            
for i in range(len(groupLabels)):

    # resample
    numberOfElements=int(1e5)
    workingDist=expressionDistributions[groupLabels[i]]
    measuredAverage=numpy.mean(workingDist)

    averageDist=[]
    for j in range(numberOfElements):
        sample=numpy.random.choice(expressionDistributions['all'],len(workingDist))
        average=numpy.mean(sample)
        averageDist.append(average)
    
    # hypothesis test
    higherRandoms=sum(numpy.greater(averageDist,measuredAverage))
    if higherRandoms > numberOfElements/2:
        pvalue=1-(higherRandoms/float(numberOfElements))
    else:
        pvalue=higherRandoms/float(numberOfElements)
    print('Group label {} has a deviation whose p-value is {}. Out of {} trials'.format(groupLabels[i],pvalue,numberOfElements))
    
    # make a figure of the expected group distribution
    x,y=histogrammer(averageDist)
    matplotlib.pyplot.plot(x,y,linestyle=':',color=theColors[i],lw=1,alpha=0.5)
    # make the line
    matplotlib.pyplot.axvline(x=measuredAverage,color=theColors[i],linestyle='-',lw=2)
    # make the dist of group

matplotlib.pyplot.xlim([-0.1,4.])
matplotlib.pyplot.ylim([-0.01,0.6])

matplotlib.pyplot.xlabel('mRNA (log$_{10}$ TPM+1)')
matplotlib.pyplot.ylabel('Probability')

matplotlib.pyplot.tight_layout()

#matplotlib.pyplot.savefig('figure.expression.distribution.TL.pdf')
#matplotlib.pyplot.savefig('figure.expression.distribution.TC.pdf')
matplotlib.pyplot.savefig('figure.expression.distribution.yellow.pdf')


matplotlib.pyplot.clf()
