###
### this script tests the hypothesis that protein expression can be explained by trends in RPF, i.e., TLR
###

import os,sys,numpy
import matplotlib,matplotlib.pyplot

#from lmfit import Model, CompositeModel

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def individualPlotter(TP,PO,TO,RO,deltaRA,RHO,TLR):

    markers=['o','s','^']
    theAlpha=0.25

    theColor='blue'
    matplotlib.pyplot.plot(TP,numpy.mean(PO,1),'-',color=theColor,label='pt',lw=0.5)
    for j in range(PO.shape[1]):
        matplotlib.pyplot.plot(TP,PO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='red'
    matplotlib.pyplot.plot(TP,numpy.mean(TO,1),'-',color=theColor,label='mRNA',lw=0.5)
    for j in range(TO.shape[1]):
        matplotlib.pyplot.plot(TP,TO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='orange'
    matplotlib.pyplot.plot(TP,numpy.mean(RO,1),'-',color=theColor,label='RBF',lw=0.5)
    for j in range(RO.shape[1]):
        matplotlib.pyplot.plot(TP,RO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='green'
    matplotlib.pyplot.plot(TP,deltaRA,'-',color=theColor,label='$\Delta$RA',lw=1)

    theColor='purple'
    matplotlib.pyplot.plot(TP,RHO,'-',color=theColor,label='$\\rho$',lw=1)

    theColor='black'
    matplotlib.pyplot.plot(TP,TLR,'-',color=theColor,label='TLR',lw=2)

    if name in ribosomalProteinGenes:
        figureName='figures/ribosomalProteins/{}.pdf'.format(name)
    else:
        figureName='figures/regularProteins/{}.pdf'.format(name)

    matplotlib.pyplot.xlabel('time (h)')
    matplotlib.pyplot.ylabel('log$_2$ FC')
    matplotlib.pyplot.title(name)
    matplotlib.pyplot.legend(loc=2)
    matplotlib.pyplot.xlim([-10,45])
    
    allValues=[PO,TO,RO,deltaRA,RHO,TLR]
    topValues=[numpy.amax(element) for element in allValues]
    bottomValues=[numpy.amin(element) for element in allValues]
    if max(topValues) > 1.9:
        a=max(topValues)+0.1*max(topValues)
    else:
        a=2
    if min(bottomValues) < -1.9:
        b=min(bottomValues)+0.1*min(bottomValues)
    else:
        b=-2    
    matplotlib.pyplot.ylim([b,a])
    
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

def proteomicsReader():

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC
    significance={} # same as data
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in data.keys():
            data[condition]={}; significance[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}; significance[condition][replicate]={}
            
        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}; significance[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                geneName=vector[0]

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                d=float(vector[2+2])
                e=float(vector[6+2])
                f=float(vector[10+2])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c

                significance[condition][replicate]['tp2vs1'][geneName]=d
                significance[condition][replicate]['tp3vs1'][geneName]=e
                significance[condition][replicate]['tp4vs1'][geneName]=f
            
    return data,significance

def riboListReader():

    riboProteins=[]
    with open(ribosomalProteinsList,'r') as f:
        for line in f:
            vector=line.split(';')
            name=vector[1].split('=')[1]
            formattedName=name.replace('_','')
            riboProteins.append(formattedName)
    uniqueList=list(set(riboProteins))

    return uniqueList

def staticAnalysis_RNAprotein():

    print('analyzing mRNA vs protein relationship in a static manner...')

    for ptReplicate in log2proteome['lysate'].keys():
        for ptTimepoint in log2proteome['lysate'][ptReplicate].keys():

            matplotlib.pyplot.plot([-7,7],[-7,7],ls='--',lw=2,color='blue')
            matplotlib.pyplot.plot([0,0],[-7,7],ls='--',lw=2,color='blue')

            for name in consistentNames:
                if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys():

                    ptRatio=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                    ptSignificance=proteomeSignificance['lysate'][ptReplicate][ptTimepoint][name]
                    mRNAratio=log2transcriptome['trna'][ptReplicate][ptTimepoint][name]

                    if ptSignificance > 0.05:
                        theColor='black'

                    else:
                        theColor='red'
                        matplotlib.pyplot.plot(mRNAratio,ptRatio,'o',alpha=0.2,mew=0,color=theColor)

                    #matplotlib.pyplot.plot(mRNAratio,ptRatio,'o',alpha=0.2,mew=0,color=theColor)

            matplotlib.pyplot.xlim([-8,8])
            matplotlib.pyplot.ylim([-8,8])

            matplotlib.pyplot.xlabel('log$_2$ FC mRNA')
            matplotlib.pyplot.ylabel('log$_2$ FC protein')

            matplotlib.pyplot.tight_layout()

            matplotlib.pyplot.savefig('figures/relation.{}.{}.png'.format(ptReplicate,ptTimepoint))
            matplotlib.pyplot.clf()

    return None

def transcriptomeRelativeConverter():

    '''
    this function assumes consistent gene names for all conditions
    '''

    log2transcriptome={}
    
    for fraction in rnaExpression.keys():
        for replicate in rnaExpression[fraction].keys():
            for timepoint in rnaExpression[fraction][replicate].keys():
                if timepoint != 'tp.1':
                    for name in rnaExpression[fraction][replicate][timepoint].keys():
                        a=rnaExpression[fraction][replicate][timepoint][name]
                        b=rnaExpression[fraction][replicate]['tp.1'][name]

                        fc=(a+1)/(b+1)
                        log2FC=numpy.log2(fc)

                        stripped=timepoint.replace('.','')
                        relativeTimePoint=stripped+'vs1'

                        if fraction not in log2transcriptome.keys():
                            log2transcriptome[fraction]={}
                        if replicate not in log2transcriptome[fraction].keys():
                            log2transcriptome[fraction][replicate]={}
                        if relativeTimePoint not in log2transcriptome[fraction][replicate].keys():
                            log2transcriptome[fraction][replicate][relativeTimePoint]={}

                        log2transcriptome[fraction][replicate][relativeTimePoint][name]=log2FC

    return log2transcriptome

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

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
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data

# 0. user defined variables

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
ribosomalProteinsList='/Volumes/omics4tb/alomana/projects/TLR/data/genome/ribosomalProteins.txt'
#transcriptomicsDataFile='/Users/adriandelomana/tmp/data/expression/expressionMatrix.kallisto.txt'
#proteomicsDataFolder='/Users/adriandelomana/tmp/data/proteomics/all/'

# 0.2. variables
timepoints=[14.3,21.5,28.8,40.8] # needs to be double checked with arjun
sortedReplicateLabels=['br1','br2','br3']
sortedRelativeTimePointLabels=['tp2vs1','tp3vs1','tp4vs1']

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression=transcriptomicsReader()
log2transcriptome=transcriptomeRelativeConverter()

# 1.2. reading protein data
log2proteome,proteomeSignificance=proteomicsReader()

# 1.3. checking consistency of transcriptome and proteome names
transcriptomeNames=[]
for fraction in log2transcriptome.keys():
    for replicate in log2transcriptome[fraction].keys():
        for timepoint in log2transcriptome[fraction][replicate].keys():
            for name in log2transcriptome[fraction][replicate][timepoint].keys():
                if name not in transcriptomeNames:
                    transcriptomeNames.append(name)

proteomeNames=[]
for fraction in log2proteome.keys():
    for replicate in log2proteome[fraction].keys():
        for timepoint in log2proteome[fraction][replicate].keys():
            for name in log2proteome[fraction][replicate][timepoint].keys():
                if name not in proteomeNames:
                    proteomeNames.append(name)

print('\t found expression quantification for {} transcripts.'.format(len(transcriptomeNames)))
print('\t found expression quantification for {} proteins.'.format(len(proteomeNames)))

# define which proteins do not have transcript equivalent
consistentNames=[]
inconsistentNames=[]
for ptName in proteomeNames:
    if ptName in transcriptomeNames:
        consistentNames.append(ptName)
    else:
        inconsistentNames.append(ptName)
consistentNames.sort()
print('\t found transcriptome info for {} proteins.'.format(len(consistentNames)))
print('\t lost {} proteins for annotation discrepancies:'.format(len(inconsistentNames)))
print('\t\t {}'.format(','.join(inconsistentNames)))
print()

# 1.4. reading ribosomal proteins list
ribosomalProteinGenes=riboListReader()

# 2. building a figure of log2 mRNA versus log2 pt
#staticAnalysis_RNAprotein()

# check how many pt are good replicates, same for transcripts. Plot it as a % of detected. use scipy.stats.norm.interval(2/3) as a rule for good replicates.

# 3. computing TLR
print('computing variables for TLR...')
TP=numpy.array(timepoints)

# 3.1. computing pt vs mRNa ratio (PMR)
rankFullProteinTrajectories=0

for name in consistentNames:

    # defining protein values
    proteinObservations=[]
    for ptReplicate in sortedReplicateLabels:
        series=[0] # manually adding the first value, log2 FC = 0 for first timepoint
        for ptTimepoint in sortedRelativeTimePointLabels:
            if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys(): 
                value=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                series.append(value)
        proteinObservations.append(series)
    PO=numpy.array(proteinObservations).T

    if PO.shape == (4,3): # 4 timepoints, 3 trajectories, full data set
        rankFullProteinTrajectories=rankFullProteinTrajectories+1

        # defining transcript values
        rnaObservations=[]
        for replicate in sortedReplicateLabels:
            series=[0]
            for tp in sortedRelativeTimePointLabels:
                value=log2transcriptome['trna'][replicate][tp][name]
                series.append(value)
            rnaObservations.append(series)
        TO=numpy.array(rnaObservations).T

        # defining RBF values
        rbfObservations=[]
        for replicate in sortedReplicateLabels:
            series=[0]
            for tp in sortedRelativeTimePointLabels:
                value=log2transcriptome['rbf'][replicate][tp][name]
                series.append(value)
            rbfObservations.append(series)
        RO=numpy.array(rbfObservations).T

        # fitting data to a model
        # https://lmfit.github.io/lmfit-py/model.html
        #mod  = CompositeModel(Model(jump), Model(gaussian), convolve)

        # computing TLR
        deltaRA=numpy.mean(RO,1)-numpy.mean(TO,1)
        RHO=numpy.mean(PO,1)-numpy.mean(TO,1)
        TLR=RHO-deltaRA
        
        # plotting
        #individualPlotter(TP,PO,TO,RO,deltaRA,RHO,TLR)

        #sys.exit()
print(rankFullProteinTrajectories)

# 3.2. testing for distribution of affinity in first points
x=[]
y=[]

for name in transcriptomeNames:
    transcriptValues=[]
    footprintValues=[]
    for replicate in sortedReplicateLabels:
        transcriptValues.append(numpy.log2(rnaExpression['trna'][replicate]['tp.1'][name]+1))
        footprintValues.append(numpy.log2(rnaExpression['rbf'][replicate]['tp.1'][name]+1))

    mT=numpy.median(transcriptValues)
    mF=numpy.median(footprintValues)
    x.append(numpy.median(transcriptValues))
    y.append(numpy.median(footprintValues))

D=numpy.array([x,y])

matplotlib.pyplot.scatter(D)
matplotlib.pyplot.xlabel('mRNA log2 (TPM+1)')
matplotlib.pyplot.ylabel('RBF log2 (TPM+1)')
matplotlib.pyplot.savefig('figure.pdf')
                
        
