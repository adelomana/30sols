###
### this script tests the hypothesis that protein expression can be explained by trends in RPF, i.e., TLR
###

import os,sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def proteomicsReader():

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

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

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c
            
    return data

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
#transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'
#proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'

transcriptomicsDataFile='/Users/adriandelomana/tmp/data/expression/expressionMatrix.kallisto.txt'
proteomicsDataFolder='/Users/adriandelomana/tmp/data/proteomics/all/'

#scipy.stats.norm.interval(2/3)

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

print('found expression quantification for {} transcripts.'.format(len(transcriptomeNames)))
print('found expression quantification for {} proteins.'.format(len(proteomeNames)))

# define which proteins do not have transcript equivalent
consistentNames=[]
inconsistentNames=[]
for ptName in proteomeNames:
    if ptName in transcriptomeNames:
        consistentNames.append(ptName)
    else:
        inconsistentNames.append(ptName)

print('found transcriptome info for {} proteins.'.format(len(consistentNames)))
print('inconsistent protein annotation for {} proteins:'.format(len(inconsistentNames)))
print(inconsistentNames)
print()

# 2. building a figure of log2 mRNA versus log2 pt
print('analyzing mRNA vs protein relationship...')

for ptReplicate in log2proteome['lysate'].keys():
    for ptTimepoint in log2proteome['lysate'][ptReplicate].keys():
        x=[]
        y=[]
        for name in consistentNames:
            if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys():
                
                ptRatio=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                mRNAratio=log2transcriptome['trna'][ptReplicate][ptTimepoint][name]
                
                x.append(mRNAratio)
                y.append(ptRatio)

        print(ptReplicate,ptTimepoint)
        print(min(x),max(x))
        print(min(y),max(y))
        print()
        matplotlib.pyplot.plot(x,y,'o',alpha=0.1,mew=0,color='black')
        

        matplotlib.pyplot.plot([-7,7],[-7,7],ls='--',lw=2,color='red')
        matplotlib.pyplot.plot([0,0],[-7,7],ls='--',lw=2,color='red')
        
        matplotlib.pyplot.xlim([-8,8])
        matplotlib.pyplot.ylim([-8,8])

        matplotlib.pyplot.xlabel('log$_2$ FC mRNA')
        matplotlib.pyplot.ylabel('log$_2$ FC protein')

        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.tight_layout()

        matplotlib.pyplot.savefig('figures/relation.{}.{}.png'.format(ptReplicate,ptTimepoint))
        matplotlib.pyplot.clf()
print()               
                

# following the steps

# 3. computing TLR
