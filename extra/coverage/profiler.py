###
### This script retrieves the coverage profiles of RNA-seq and Ribo-seq for all ribosomal protein genes. It stores it as text files.
###

import sys,numpy,HTSeq
import multiprocessing,multiprocessing.pool

def analysis(genomicFeature):

    '''
    This function computes the histograms of reads across transcript lengths.
    '''

    print('\t computing coverage for {}...'.format(genomicFeature))

    
    # f.1 define window of coverage depending if it's an operon or a gene
    print('\t\t computing window...')
    if genomicFeature in riboOperons.keys(): # work with operons

        print(genomicFeature)

        # obtain the relevant features
        contigs=[]; starts=[]; ends=[]; strands=[]
        localGenes=riboOperons[genomicFeature]

        for feature in annotationObject:
            if feature.type == 'gene':
                strippedID=feature.attr['ID'].replace('_','')
                if strippedID in localGenes:
                    
                    contig=feature.iv.chrom
                    start=feature.iv.start-margin
                    end=feature.iv.end+margin
                    strand=feature.iv.strand

                    contigs.append(contig); starts.append(start); ends.append(end); strands.append(strand)

        # check consistency of strands
        if len(list(set(strands))) > 1:
            print('Detected gene in operon with different orientation. Exiting...')
            sys.exit()
            
        # define positions for coverage computing
        contig=contigs[0]
        start=min(starts)-margin
        end=max(ends)+margin
        strand=strands[0]

        print(start,end)

        windowP=HTSeq.GenomicInterval(contig,start,end,"+")
        windowM=HTSeq.GenomicInterval(contig,start,end,"-")

    else: # work with genes
        for feature in annotationObject:
            if feature.type == 'gene':
                strippedID=feature.attr['ID'].replace('_','')
                if strippedID == genomicFeature:
                    break
        # define positions for coverage computing
        print(feature)
        contig=feature.iv.chrom
        start=feature.iv.start-margin
        end=feature.iv.end+margin
        strand=feature.iv.strand

        windowP=HTSeq.GenomicInterval(contig,start,end,"+")
        windowM=HTSeq.GenomicInterval(contig,start,end,"-")

    """
    # f.2. compute coverage based on window
    print('\t\t computing coverage...')
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
                fileName='{}{}.{}.{}.{}.txt'.format(coverageDir,timepoint,replicate,genomicFeature,experiment)
                f=open(fileName,'w')
                f.write('# name {}\n'.format(genomicFeature))
                f.write('# timepoint {}\n'.format(timepoint))
                f.write('# replicate {}\n'.format(replicate))
                f.write('# strand {}\n'.format(strand))
                f.write('# experiment {}\n'.format(experiment))
                f.write('# sumP,sumM {},{}\n'.format(sum(profileP),sum(profileM)))
                f.write('# location \t counts on strand plus \t counts on strand minus\n')
                for i in range(len(pos)):
                    f.write('{}\t{}\t{}\n'.format(pos[i],profileP[i],profileM[i]))
                f.close()
    """

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
bamFilesDir='/proj/omics4tb/alomana/projects/TLR/data/BAM/'
annotationFile='/proj/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'
coverageDir='/proj/omics4tb/alomana/projects/TLR/data/coverage/'
operonPredictionsDir='/proj/omics4tb/alomana/projects/TLR/data/microbesOnline/'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

margin=100 # excess of base pairs

# 1. read data
print('Reading data...')
riboOperons,NORPGs=dataReader()

# 2. iterate analysis over ribosomal proteins
print('Performing analysis...')

# 2.1. read annotation file
annotationObject=HTSeq.GFF_Reader(annotationFile)

# 2.2. selecting appropriate genomic locations
genomicFeatures=list(riboOperons.keys())+NORPGs

# 2.3.a. iterate over genomicFeatures in a parallel manner
numberOfThreads=len(genomicFeatures)
print('Initialized parallel analysis using {} threads...'.format(numberOfThreads))
hydra=multiprocessing.pool.Pool(numberOfThreads)
tempo=hydra.map(analysis,genomicFeatures)
print('... completed.')

# 2.3.b. iterate over genomicFeatures single-thread
for genomicFeature in genomicFeatures:
#for genomicFeature in NORPGs:
    analysis(genomicFeature)
