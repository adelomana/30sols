###
### This script looks for different GREs between two sets of transcripts. need not to bump into another gene.
###

import sys
import matplotlib,matplotlib.pyplot

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
        ###limits=[geneHead-100,geneHead]
                
    elif geneStrand == '-': # what is the first 3' head?
        allHeads=[]
        for element in allGeneCoordinates:
            head=allGeneCoordinates[element][0]
            if head > geneTail:
                allHeads.append(head)
        allHeads.sort()
        curb=allHeads[0]
        limits=[geneTail,curb]
        ###limits=[geneTail,geneTail+100]

    else:
        print('error when selecting strand')
        sys.exit()
                    
    # f.2. cut the right sequence
    room=limits[1]-limits[0]

    print(limits,room)
    
    fullUpstream=genomeSequence[limits[0]-1:limits[1]-1]

    # f.3. trim sequence
    if geneStrand == '+':
        if len(fullUpstream) <= upstreamSearch:
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
groupingDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/ribo.groupings.csv'
riboOperonsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

# 1. read data

# 1.1. read info about two groups of transcripts
geneSets={}
geneSets['group A']=[]
geneSets['group B']=[]
expressionCoordinates={}

with open(groupingDataFile,'r') as f:
    next(f)
    for line in f:
        v=line.split(',')
        if v[0] == 'group A':
            geneSets[v[0]].append(v[1])
        else:
            geneSets[v[0]].append(v[1])

        expressionCoordinates[v[1]]=[float(v[2]),float(v[3])]

# 1.2. read operons
rbptOperons={}
genesInOperons=[]
with open(riboOperonsFile,'r') as f:
    next(f)
    for line in f:
        v=line.split()
        name=v[0]
        rbptOperons[name]=[]
        
        elements=v[1:]
        for element in elements:
            rbptOperons[name].append(element)
            genesInOperons.append(element)

# 1.3. read gene orientations
geneOrientations={}
with open(annotationFile,'r') as f:
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        if 'gene' in v:
            info=v[-1]

            strand=v[6]
            id=v[8].split('ID=')[1].split(';')[0]
            geneOrientations[id]=strand

# 1.4 convert new annotation to old annotation
annotationMap={}
with open(annotationFile,'r') as f:
    next(f)
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        print(v)
        sys.exit()
sys.exit()

# 2. trim gene sets to operon heads
geneLeaders=[] # genes that are either not part of an operon or operon headers

# 2.1. check if gene is in operon
for element in expressionCoordinates:
    if element not in genesInOperons:
        geneLeaders.append(element)
    else:
        # 2.2. if gene is in an operon, make sure it is at its head
        if geneOrientations[element] == '+':
            for operon in rbptOperons:
                if element == rbptOperons[operon][0]:
                    geneLeaders.append(element)
        if geneOrientations[element] == '-':
            for operon in rbptOperons:
                if element == rbptOperons[operon][-1]:
                    geneLeaders.append(element)


print(rbptOperons)
print()
print(geneLeaders)

# 3. plot fold-changes of gene leaders
for leader in geneLeaders:
    x=expressionCoordinates[leader][1]
    y=expressionCoordinates[leader][0]
    matplotlib.pyplot.plot(x,y,'ok')

matplotlib.pyplot.xlabel('median expression')
matplotlib.pyplot.ylabel('fold-change')
matplotlib.pyplot.savefig('figure.gene.expression.leaders.pdf')

# 4. obtain GRE counts for each gene

# 4.1. convert from new annotation to old annotation

# 4.2. define the region of interest considering the gff3 file provided by Wei-ju

# 4.3. convert each GRE PDF into an integral


# 5. make a PCA plot of proportions, or face values of integrals



sys.exit()
                    
           
# 2. define the upstream regulatory sequence for each gene, independently of being inside an operon
upstreamSections={}
upstreamSections['groupA']={}
upstreamSections['groupB']={}

regulatoryReferences={}

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
        if referenceGene not in upstreamSections[geneSetName]:
            upstreamSections[geneSetName][referenceGene]=[upstreamRegion,[geneName]]
        else:
            upstreamSections[geneSetName][referenceGene][1].append(geneName)
        
for geneSetName in upstreamSections:
    for reference in upstreamSections[geneSetName]:
        
        print('>{}.{}.{}.{}'.format(geneSetName,reference,'-'.join(upstreamSections[geneSetName][reference][1]),geneCoordinates[reference][2]))
        print(upstreamSections[geneSetName][reference][0])

    print()
        

# 4. detect GREs using MEME
