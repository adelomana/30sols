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
    pos=[]; cov=[]
    with open(dataFileName,'r') as f:
        for line in f:
            vector=line.split()

            # f.1. obtaining the metadata
            if vector[1] == 'strand':
                strand=vector[2]
            if vector[1] == 'experiment':
                experiment=vector[2]

            #if vector[1] == 'sumP,sumM':
            #    print(vector)

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
                cov.append(int(vector[column]))

    # dealing with positions
    if strand == '-':
        pos=pos[::-1]
    p=numpy.array(pos)
    normalizedPosition=p-min(p)

    print(name,synonyms[name],strand)

    return normalizedPosition,cov

def figureMaker():

    '''
    This function builds a figure of the coverage of reads over ribo-pt genes.
    '''

    # f.1. read the data
    for experiment in experiments:
        for timepoint in timepoints:
            y=[]
            for replicate in replicates:
                dataFileName='{}{}.{}.{}.{}.txt'.format(coverageDir,timepoint,replicate,name,experiment)
                pos,cov=coverageFileReader(dataFileName)
                y.append(cov)

            

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
        
        matplotlib.pyplot.title('{} {}'.format(synonyms[name],flag))

        #matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=1,ncol=2,fontsize=14)

        figureName='figures/figure.{}.{}.pdf'.format(synonyms[name],experiment)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

    return None

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    synonyms={}
    
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])
            synonyms[vector[0]]=vector[2].replace('\n','')

    riboPtNames.sort()
            
    return riboPtNames,synonyms

###
### MAIN
###

# 0. user defined variables
coverageDir='/Volumes/omics4tb/alomana/projects/TLR/data/coverage/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

colors=['red','orange','green','blue']

# 1. iterate over ribo-pt genes
riboPtNames,synonyms=riboPtNamesReader()

# 2. build figure
for name in riboPtNames:
    print('building figure for {}...'.format(name))
    figureMaker()
