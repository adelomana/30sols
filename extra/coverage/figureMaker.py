###
### This script builds histograms from the coverage profile text files.
###

import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def coverageFileReader(dataFileName):

    '''
    This function reads the coverage files and returns the positions and the coverage.
    '''

    strand=None; experiment=None
    pos=[]; coverage=[]
    with open(dataFileName,'r') as f:
        for line in f:
            vector=line.split()

            # f.1. obtaining the metadata
            if vector[1] == 'strand':
                strand=vector[2]
            if vector[1] == 'experiment':
                experiment=vector[2]

            # f.2. reading the information
            if vector[0] != '#':
                
                # define which column to read
                if experiment == 'rbf':
                    if strand == '+':
                        column=1
                    elif strand == '-':
                        column=2
                    else:
                        print('Error selecting strand at rbf. Exiting...')
                        sys.exit()
                
                elif experiment == 'trna':
                    if strand == '+':
                        column=2
                    elif strand == '-':
                        column=1
                    else:
                        print('Error selecting strand at trna. Exiting...')
                        sys.exit()
                else:
                        print(experiment)
                        print('Error from experiment selection. Exiting...')
                        sys.exit()
                        
                # read columns
                pos.append(int(vector[0]))
                coverage.append(int(vector[column]))

    # dealing with positions
    if strand == '-':
        pos=pos[::-1]
    p=numpy.array(pos)
    normalizedPosition=p-min(p)-margin

    print(genomicFeature,strand)

    return normalizedPosition,coverage

def figureMaker(genomicFeature):

    '''
    This function builds a figure of the coverage of reads over genomic features.
    '''

    # f.1. read the data
    for experiment in experiments:
        for timepoint in timepoints:
            y=[]
            for replicate in replicates:
                dataFileName='{}{}.{}.{}.{}.txt'.format(coverageDir,timepoint,replicate,genomicFeature,experiment)
                pos,coverage=coverageFileReader(dataFileName)
                y.append(coverage)

            # compute PDF 
            average=numpy.mean(numpy.array(y),axis=0)
            pdf=average/sum(average)

            # define the color
            theColor=colors[timepoints.index(timepoint)]

            # plot
            matplotlib.pyplot.plot(pos,pdf,'-',color=theColor,label=timepoint)

        # f.2.3 final figure closing
        matplotlib.pyplot.xlabel("Relative genomic position (5'->3')")
        matplotlib.pyplot.ylabel('p(coverage)')
        if experiment == 'trna':
            flag='RNA-seq'
        else:
            flag='Ribo-seq'
        
        matplotlib.pyplot.title('{} {}'.format(genomicFeature,flag))

        #matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=1,ncol=2,fontsize=14)

        figureName='figures/figure.{}.{}.pdf'.format(genomicFeature,experiment)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

    return None

def dataReader():

    '''
    This function reads the ribosomal protein operons and genes.
    '''

    # f.1. ribo-pt gene operons
    operonPredictions={}
    fileName=operonPredictionsDir+'riboPtOperons.txt'
    with open(fileName,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            name=vector[0]
            genes=[]
            for i in range(len(vector)-1):
                gene=vector[i+1].replace('\n','')
                genes.append(gene)
            operonPredictions[name]=genes

    # f.2. non-operon ribo-pt genes
    NORPGs=[]
    fileName=operonPredictionsDir+'NORPGs.txt'
    with open(fileName,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            name=vector[0].replace('\n','')
            NORPGs.append(name)

    # f.3. print information about retrieval
    a=[]
    for operon in operonPredictions:
        for name in operonPredictions[operon]:
            if name not in a:
                a.append(name)
    print('\t Recovered {} genes in {} operons.'.format(len(a),len(operonPredictions)))
    print('\t Recovered {} genes not in operons.'.format(len(NORPGs)))
    for name in NORPGs:
        if name not in a:
            a.append(name)
    print('\t Total genes recovered: {}'.format(len(a)))
            
    return operonPredictions,NORPGs

###
### MAIN
###

# 0. user defined variables
coverageDir='/Volumes/omics4tb/alomana/projects/TLR/data/coverage/'
operonPredictionsDir='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

colors=['red','orange','green','blue']

margin=100 # excess of base pairs

# 1. read data
print('Reading data...')
riboOperons,NORPGs=dataReader()
genomicFeatures=list(riboOperons.keys())+NORPGs

# 2. build figure
for genomicFeature in genomicFeatures:
    print('building figure for {}...'.format(genomicFeature))
    figureMaker(genomicFeature)
print('... completed.')
