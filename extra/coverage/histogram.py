###
### This script builds histograms of RNA-seq and Ribo-seq for all ribosomal protein genes
###

import sys,numpy,HTSeq
import matplotlib,matplotlib.pyplot
import multiprocessing,multiprocessing.pool

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def analysis():

    '''
    This function computes the histograms of reads across transcript lengths.
    '''




    #rbf.rep.1.tp.1/Aligned.sortedByCoord.out.bam
    
    localFeature=genomicFeatures[riboPt]
    
    print(localFeature)

    contig=localFeature.iv.chrom
    start=localFeature.iv.start
    end=localFeature.iv.end
    strand=localFeature.iv.strand

    print(start)
    print(end)
    print(contig)
    print(strand)

    window=HTSeq.GenomicInterval(contig,start,end,".")

    coverage=HTSeq.GenomicArray("auto",stranded=False,typecode="i")
    

    
    for timepoint in timepoints:
        for replicate in replicates:

            # f.1. define the bam file
            bamFile=bamFilesDir+'rbf.{}.{}/Aligned.sortedByCoord.out.bam'.format(replicate,timepoint)

            # f.2. read BAM file
            sortedBAMfile=HTSeq.BAM_Reader(bamFile)
            for alignment in sortedBAMfile:
                if alignment.aligned:
                    #alignment.iv.length = 50 ###### double check if this makes any difference!!!!
                    coverage[ alignment.iv ] += 1

            profile=list(coverage[window])

            print(profile,sum(profile))


            
            loc=numpy.arange(start,end)
            print('loc',loc)
            if strand == '+':
                pos=loc
            elif strand == '-':
                pos=loc[::-1]
            else:
                print('error at strand selection')
                sys.exit()

            print('pos',pos)

            # writing a file


            # making a plot
            matplotlib.pyplot.plot(pos,profile,'-b')


            figureFile='figure.pdf'

            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig(figureFile)
            matplotlib.pyplot.clf()


            ##### inverse
            matplotlib.pyplot.plot(pos[::-1],profile,'-r')
            matplotlib.pyplot.tight_layout()
            
            matplotlib.pyplot.savefig('inverse.pdf')
            matplotlib.pyplot.clf()
            sys.exit()

            


    # compute the pdf over transcript length

    # compute the normalized amount (with respect to library size) of reads across transcript length

    return None

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

###
### MAIN
###

# 0. user defined variables
bamFilesDir='/Volumes/omics4tb/alomana/projects/TLR/data/BAM/'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']

# 1. read data
print('reading data...')
riboPtNames=riboPtNamesReader()

# 2. iterate analysis over ribosomal proteins
print('performing analysis...')

# 2.1. read annotation file
annotationObject=HTSeq.GFF_Reader(annotationFile)

# 2.2. selecting appropriate genomic locations
genomicFeatures={}
for feature in annotationObject:
    if feature.type == 'gene':
        strippedID=feature.attr['ID'].replace('_','')
        if strippedID in riboPtNames:
            genomicFeatures[strippedID]=feature

# 2.2. iterate over genes

##### this section needs to be paralelized
for riboPt in riboPtNames:
    print('building figure for {}...'.format(riboPt))
    analysis()



