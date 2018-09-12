###
### This script looks for different GREs between two sets of transcripts. need not to bump into another gene.
###

import sys

def fastaFileReader():

    '''
    This function reads the fasta file for the genome and returns a large string.
    '''

    genomeSequence=''

    completed=False
    with open(genomeFile,'r') as f:
        next(f)
        for line in f:
            v=line.split()
            if len(v) != 0:
                genomeSequence=genomeSequence+v[0]
            else:
                break
    
    return genomeSequence

def upstreamSelector():

    '''
    This function appends to "geneCoordinates" variable the upstream region.
    '''    

    for geneSetName in geneSets:
        for geneName in geneSets[geneSetName]:

            

            # f.1. find the sequence up to the next gene
            geneHead=geneCoordinates[geneName][0]
            geneTail=geneCoordinates[geneName][1]
            candidateStrand=geneCoordinates[geneName][2]

            print('working with',geneSetName,geneName,candidateStrand)
       
            if candidateStrand == '+': # what is the gene tail smaller than the gene head?
                allTails=[]
                for element in allGeneCoordinates:
                    tail=allGeneCoordinates[element][1]
                    if tail < geneHead:
                        allTails.append(tail)
                allTails.sort()
                curb=allTails[-1]
                limits=[curb,geneHead]
                
            else: # what is the first 3' head?
                allHeads=[]
                for element in allGeneCoordinates:
                    head=allGeneCoordinates[element][0]
                    if head > geneTail:
                        allHeads.append(head)
                allHeads.sort()
                curb=allHeads[0]
                limits=[geneTail,curb]
                
            print(geneName,limits)
            print()
                
        


        # f.2. cut the right sequence
        #if len(fullUpstream) > upstreamWindow:
        #    specificUpstream=fullUpstream[:upstreamWindow]
        #else:
        #    specificUpstream=fullUpstream

    return geneCoordinates

### MAIN

# 0. user defined variables

# 0.1. paths
groupingDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/data.txt'
genomeFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.fasta'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
riboOperonsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'

# 0.2. options
upstreamWindow=50

# 1. read data

# 1.1. read info about two groups of transcripts
geneSets={}
geneSets['groupA']=[]
geneSets['groupB']=[]
allGenes=[]

with open(groupingDataFile,'r') as f:
    next(f)
    for line in f:
        v=line.split('\t')
        groupTag=v[0].replace(' ','')
        geneName=v[1].replace(' ', '')
        geneSets[groupTag].append(geneName)
        allGenes.append(geneName)
        
# 1.2. read genome info
geneCoordinates={}
recoveredGenes=[]
allGeneCoordinates={}
synonyms={}

with open(annotationFile,'r') as f:
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        if 'gene' in v:
            info=v[-1]
            start=int(v[3])
            stop=int(v[4])
            id=v[8].split('ID=')[1].split(';')[0]
            allGeneCoordinates[id]=[start,stop]

            if 'old_locus_tag=' in info:
                geneName=info.split('old_locus_tag=')[1].split('%')[0].replace('\n','')
                if geneName in allGenes:
                    recoveredGenes.append(geneName)

                    contig=v[0]
                    strand=v[6]
                    geneCoordinates[geneName]=[start,stop,strand]

                    synonyms[geneName]=id

# check about full recovery of elements
if len(recoveredGenes) != len(allGenes):
    print('error at recovering genes...')
    print(len(recoveredGenes),len(allGenes))
    sys.exit()

# 1.3. read FASTA file
genomeSequence=fastaFileReader()

# 2. manipulate data: remove genes that are behind an operon head
operonFollowers=[]
with open(riboOperonsFile,'r') as f:
    next(f)
    for line in f:
        v=line.split()
        followers=v[2:]
        for follower in followers:
            operonFollowers.append(follower)
print(len(operonFollowers),'operon followers found.')

# iterate over gene sets: check if they are operon followers
cursed=[] 
for geneSetName in geneSets:
    for geneName in geneSets[geneSetName]:
        candidate=synonyms[geneName]
        if candidate in operonFollowers:
            geneSets[geneSetName].remove(geneName)

# 3. define for each gene the regulatory region
regulatorySequences=upstreamSelector()


# 4. detect GREs using MEME
