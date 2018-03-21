###
### This script retrieves the coverage profiles of RNA-seq and Ribo-seq for all ribosomal protein genes. It stores it as text files.
###

import sys,numpy,HTSeq
import multiprocessing,multiprocessing.pool

def analysis(riboPt):

    '''
    This function computes the histograms of reads across transcript lengths.
    '''

    print('\t building figure for {}...'.format(riboPt))
    
    localFeature=genomicFeatures[riboPt]
    
    contig=localFeature.iv.chrom
    start=localFeature.iv.start
    end=localFeature.iv.end
    strand=localFeature.iv.strand

    windowP=HTSeq.GenomicInterval(contig,start,end,"+")
    windowM=HTSeq.GenomicInterval(contig,start,end,"-")

    coverage=HTSeq.GenomicArray("auto",stranded=True,typecode="i")
        
    for timepoint in timepoints:
        for replicate in replicates:
            for experiment in experiments:

                # f.1. define the bam file
                bamFile=bamFilesDir+'{}.{}.{}/Aligned.sortedByCoord.out.bam'.format(experiment,replicate,timepoint)

                # f.2. read BAM file
                sortedBAMfile=HTSeq.BAM_Reader(bamFile)
                for alignment in sortedBAMfile:
                    if alignment.aligned:
                        coverage[ alignment.iv ] += 1
                        
                # f.3. compute coverage
                profileP=list(coverage[windowP])
                profileM=list(coverage[windowM])

                # f.4. define genomic positions with respect to strands
                loc=numpy.arange(start,end)
                if strand == '+':
                    pos=loc
                elif strand == '-':
                    pos=loc[::-1]
                else:
                    print('error at strand selection')
                    sys.exit()

                # f.5. writing a file
                fileName='{}{}.{}.{}.{}.txt'.format(coverageDir,timepoint,replicate,riboPt,experiment)
                f=open(fileName,'w')
                f.write('# name {}\n'.format(riboPt))
                f.write('# timepoint {}\n'.format(timepoint))
                f.write('# replicate {}\n'.format(replicate))
                f.write('# strand {}\n'.format(strand))
                f.write('# experiment {}\n'.format(experiment))
                f.write('# sumP,sumM {},{}\n'.format(sum(profileP),sum(profileM)))
                f.write('# location \t counts on strand plus \t counts on strand minus\n')
                for i in range(len(pos)):
                    f.write('{}\t{}\t{}\n'.format(pos[i],profileP[i],profileM[i]))
                f.close()

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

    riboPtNames.sort()
            
    return riboPtNames

###
### MAIN
###

# 0. user defined variables
bamFilesDir='/proj/omics4tb/alomana/projects/TLR/data/BAM/'
annotationFile='/proj/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'
ribosomalProteinsFile='/proj/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
coverageDir='/proj/omics4tb/alomana/projects/TLR/data/coverage/'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

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

# 2.3. iterate over genes in a parallel manner
numberOfThreads=len(riboPtNames)
hydra=multiprocessing.pool.Pool(numberOfThreads)
tempo=hydra.map(analysis,riboPtNames)
