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

    # f.1. find the sequence up to the next gene
    geneHead=geneCoordinates[referenceGene][0]
    geneTail=geneCoordinates[referenceGene][1]
    geneStrand=geneCoordinates[referenceGene][2]

       
    if geneStrand == '+': # what is the gene tail smaller than the gene head?
        allTails=[]
        for element in allGeneCoordinates:
            tail=allGeneCoordinates[element][1]
            if tail < geneHead:
                allTails.append(tail)
        allTails.sort()
        curb=allTails[-1]
        limits=[curb,geneHead]
                
    elif geneStrand == '-': # what is the first 3' head?
        allHeads=[]
        for element in allGeneCoordinates:
            head=allGeneCoordinates[element][0]
            if head > geneTail:
                allHeads.append(head)
        allHeads.sort()
        curb=allHeads[0]
        limits=[geneTail,curb]

    else:
        print('error when selecting strand')
        sys.exit()
                    
    # f.2. cut the right sequence
    room=limits[1]-limits[0]
    fullUpstream=genomeSequence[limits[0]-1:limits[1]-1]

    # f.3. trim sequence
    if geneStrand == '+':
        if len(fullUpstream) < upstreamSearch:
            trimmed=fullUpstream
        elif len(fullUpstream) > upstreamSearch:
            trimmed=fullUpstream[len(fullUpstream)-upstreamSearch:]
        else:
            print('really unlucky')
            sys.exit()
    elif geneStrand == '-':
        trimmed=fullUpstream[:upstreamSearch]
    else:
        print('error when selecting strand while trimming sequence')
        sys.exit()

    return trimmed

### MAIN

# 0. user defined variables

# 0.1. paths
groupingDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/data.txt'
genomeFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.fasta'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
riboOperonsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'

upstreamSearch=250

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
geneCoordinates={}; allGeneCoordinates={}
recoveredGenes=[]
synonyms={}; invertedSynonyms={}

with open(annotationFile,'r') as f:
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        if 'gene' in v:
            info=v[-1]
            
            start=int(v[3])
            stop=int(v[4])
            strand=v[6]
            id=v[8].split('ID=')[1].split(';')[0]

            allGeneCoordinates[id]=[start,stop,strand]

            if 'old_locus_tag=' in info:
                geneName=info.split('old_locus_tag=')[1].split('%')[0].replace('\n','')
                geneCoordinates[geneName]=[start,stop,strand]

                synonyms[geneName]=id
                invertedSynonyms[id]=geneName

# 1.3. read FASTA file
genomeSequence=fastaFileReader()

# 1.4. read operons
rbptOperons={}
with open(riboOperonsFile,'r') as f:
    next(f)
    for line in f:
        v=line.split()
        name=v[0]
        rbptOperons[name]=[]
        
        elements=v[1:]
        for element in elements:
            geneID=invertedSynonyms[element]
            rbptOperons[name].append(geneID)

# 2. define the upstream regulatory sequence for each gene, independently of being inside an operon
upstreamSections={}
upstreamSections['groupA']={}
upstreamSections['groupB']={}

# for each gene, find if it belongs to an operon
for geneSetName in geneSets:
    print(geneSetName)
    for geneName in geneSets[geneSetName]:

        # define gene strand direction
        strand=geneCoordinates[geneName][2]

        # find if it belongs to an operon
        operonMembership=None
        for operon in rbptOperons:
            if geneName in rbptOperons[operon]:
                operonMembership=operon

        if operonMembership != None: # if it belongs to an operon, select the head upstream region
            if strand == '+':
                referenceGene=rbptOperons[operonMembership][0]
            elif strand == '-':
                referenceGene=rbptOperons[operonMembership][-1]
            else:
                print('error when defining strand')
                sys.exit()
            
        else: # if it does not belong to an operon, select the adjacent upstream region
            referenceGene=geneName

        # obtain the upstream region
        upstreamRegion=upstreamSelector()

        # trim region and append to dictionary            
        upstreamSections[geneSetName][referenceGene]=upstreamRegion
        
for geneSetName in upstreamSections:
    for geneName in upstreamSections[geneSetName]:
        print('>{}.{}.{}'.format(geneName,geneSetName,geneCoordinates[geneName][2]))
        print(upstreamSections[geneSetName][geneName])
    print()
        

# 4. detect GREs using MEME
