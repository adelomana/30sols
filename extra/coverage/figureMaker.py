###
### This script builds histograms from the coverage profile text files.
###

import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def cdsBlockDefiner(genomicFeature):

    '''
    This function defines the CDS blocks for each genomic feature.
    '''

    cdsBlocks=[]
    geneIDs=[]

    if genomicFeature in riboOperons:
        geneIDs=riboOperons[genomicFeature]
    else:
        geneIDs=[genomicFeature]

    for geneID in geneIDs:

        start=geneAnnotations[geneID][0]
        end=geneAnnotations[geneID][1]
        strand=geneAnnotations[geneID][2]
        cdsBlocks.append([start,end,strand,geneID])

    return cdsBlocks

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

    return normalizedPosition,coverage

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

def figureMaker(genomicFeature,cdsBlocks):

    '''
    This function builds a figure of the coverage of reads over genomic features.
    '''
    
    # f.1. iterate over experiments
    for experiment in experiments:

        fig=matplotlib.pyplot.figure()
        ax=fig.add_subplot(111)
        heaven=0
        for timepoint in timepoints:
            y=[]
            for replicate in replicates:
                dataFileName='{}{}.{}.{}.{}.txt'.format(coverageDir,timepoint,replicate,genomicFeature,experiment)
                pos,coverage=coverageFileReader(dataFileName)
                y.append(coverage)

            # compute PDF 
            average=numpy.mean(numpy.array(y),axis=0)
            pdf=average/sum(average)
            if heaven < numpy.max(pdf):
                heaven=numpy.max(pdf)

            # define the color
            theColor=colors[timepoints.index(timepoint)]

            # plot
            ax.plot(pos,pdf,'-',color=theColor,label=timepoint)

        # f.1.2. make boxes
        strand=cdsBlocks[0][2]
        boxSize=heaven*0.1
        postTexty=-boxSize/2
        
        for block in cdsBlocks:

            geneName=block[-1]
            boxLabel=synonyms[geneName][0]

            # stand-specific options
            if strand == '+':

                ref=cdsBlocks[0][0]
                start=block[0]-ref
                end=block[1]-ref

            else:
                
                ref=cdsBlocks[-1][1]
                start=ref-block[1]
                end=ref-block[0]

            # make boxes
            geneBox=matplotlib.patches.Rectangle(xy=(start,-boxSize),width=(end-start),height=boxSize,fc='white',lw=1,ec='black')
            ax.add_patch(geneBox)

            # define text
            posTextx=start+(end-start)/2
            if end-start > 500:
                theFontSize=10
            else:
                theFontSize=7
            ax.text(posTextx,postTexty,boxLabel,color='black',horizontalalignment='center',verticalalignment='center',fontsize=theFontSize)

        # f.1.3 final figure closing
        matplotlib.pyplot.xlabel("Relative genomic position (5'->3')")
        matplotlib.pyplot.ylabel('p(coverage)')
        if experiment == 'trna':
            flag='RNA-seq'
        else:
            flag='Ribo-seq'

        matplotlib.pyplot.title('{} {}'.format(genomicFeature,flag))

        matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=0,ncol=2,fontsize=14)

        matplotlib.pyplot.ylim([-1.5*boxSize,boxSize*10.5])

        figureName='figures/figure.{}.{}.pdf'.format(genomicFeature,experiment)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

    return None

def geneAnnotationReader():

    geneAnnotations={}

    with open(gffFile,'r') as f:
        next(f)
        next(f)
        for line in f:
            vector=line.split('\t')
            if len(vector) > 3:
                if vector[2] == 'gene': # check if gene and cds match exact position for all
                    name=vector[8].split(';')[0].replace('ID=','')
                    start=int(vector[3])
                    end=int(vector[4])
                    strand=vector[6]
                    geneAnnotations[name]=[start,end,strand]

    return geneAnnotations

def synonymsReader():

    '''
    This function reads the GFF3 file and returns a dictionary with synonyms between old and new locus names.
    '''

    synonyms={}
    with open(gffFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1].replace('\n','')
                if 'old_locus_tag=' in info:
                    old=info.split('old_locus_tag=')[1].split(';')[0]
                    new=info.split('ID=')[1].split(';')[0]

                    if '%' in old:
                        olds=old.split('%2C')
                        synonyms[new]=[old for old in olds]
                    else:
                        synonyms[new]=[old]
                
    return synonyms

###
### MAIN
###

# 0. user defined variables
coverageDir='/Volumes/omics4tb/alomana/projects/TLR/data/coverage/'
operonPredictionsDir='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/'
gffFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

colors=['red','orange','green','blue']

margin=100 # excess of base pairs

# 1. read data
print('reading data...')

# 1.1. read gff3 file
geneAnnotations=geneAnnotationReader()

# 1.2. read operon memberships
riboOperons,NORPGs=dataReader()
genomicFeatures=list(riboOperons.keys())+NORPGs
genomicFeatures.sort()

# 1.3. define synonyms
synonyms=synonymsReader()

# 2. build figure
print('building figures...')

for genomicFeature in genomicFeatures:
    print('building figure for {}...'.format(genomicFeature))
    cdsBlocks=cdsBlockDefiner(genomicFeature)
    figureMaker(genomicFeature,cdsBlocks)

print('... completed.')
