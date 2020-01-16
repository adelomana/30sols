###
### This script looks for different GREs between two sets of transcripts. need not to bump into another gene.
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.decomposition

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

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
regulatorySequenceSize=100

# 0.1. paths
GREsDir='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/positions/hal_gres/'
groupingDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/ribo.groupings.csv'
riboOperonsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
positionsFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/positions/hal_genes.tsv'

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
annotationMap={}; reverseAnnotationMap={}
with open(annotationFile,'r') as f:
    next(f)
    next(f)
    next(f)
    for line in f:
        v=line.split('\t')
        info=v[-1]
        if 'old_locus_tag' in info:
            new=info.split('ID=')[1].split(';')[0]
            old=info.split('old_locus_tag=')[1].split(';')[0].split('\n')[0]
            annotationMap[new]=old
            reverseAnnotationMap[old]=new

# 1.5. associate regulatory regions to genes in the old annotation
regulatoryRegions={}
with open(positionsFile,'r') as f:
    next(f)
    for line in f:
        v=line.split('\t')
        name=v[1]
        a=int(v[3])
        strand=v[5]
        if strand == '+':
            region=[a-regulatorySequenceSize,a]
        elif strand == '-':
            region=[a,a+regulatorySequenceSize]
        else:
            print('error while reading strand')
            sys.exit()
        regulatoryRegions[name]=region
        
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

# 3. plot fold-changes of gene leaders
for leader in geneLeaders:
    x=expressionCoordinates[leader][1]
    y=expressionCoordinates[leader][0]
    if leader in geneSets['group A']:
        matplotlib.pyplot.plot(x,y,'or')
    else:
        matplotlib.pyplot.plot(x,y,'ob')
        print('bLEADER',leader)

matplotlib.pyplot.xlabel('median expression')
matplotlib.pyplot.ylabel('fold-change')

matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig('figure.gene.expression.leaders.pdf')
matplotlib.pyplot.clf()

# 4. obtain GRE counts for each gene
geneGREs={}
for leader in geneLeaders:
    
    # 4.1. convert from new annotation to old annotation
    oldAnnotation=annotationMap[leader]

    print(leader,oldAnnotation)
    
    # 4.2. define the region of interest considering the gff3 file provided by Wei-ju
    regionOfInterest=regulatoryRegions[oldAnnotation]
    regulatoryDomain=numpy.arange(regionOfInterest[0],regionOfInterest[1])

    # 4.3. convert each GRE PDF into an integral above a threshold of say 10
    path2search=GREsDir+oldAnnotation
    allFiles=os.listdir(path2search)
    GREfiles=[element for element in allFiles if element[:3] == 'hal']
    for dimension in GREfiles:

        GREID=int(dimension.split('hal_')[1].split('.tsv')[0])
        dataFile=path2search+'/'+dimension


        # recover GREs above a threshold=10
        threshold=10
        counts=[]
        with open(dataFile,'r') as f:
            next(f)
            for line in f:
                v=line.split('\t')
                pos=int(v[0])
                count=int(v[1])
                hit=pos in regulatoryDomain

                if hit == True and count >= threshold:
                    counts.append(count)

        # compute the integral
        if len(counts) > 0:
            if oldAnnotation not in geneGREs:
                geneGREs[oldAnnotation]={}
            
            integral=sum(counts)/len(counts)
            geneGREs[oldAnnotation][GREID]=integral
            print('#####',GREID,integral)
    print()

# 5. make a PCA plot of proportions, or face values of integrals

# 5.1. retrieve labels of the dimensions
dimensions=[]
GREgenes=list(geneGREs.keys())
GREgenes.sort()
                  
for gene in GREgenes:
    for dimension in geneGREs[gene]:
        if dimension not in dimensions:
            dimensions.append(dimension)
dimensions.sort()

# 5.2. compute a table with GREs 9 and 28 or others, grouped by group A or group B
GREdistribution={}
GREdistribution['group A']={}
GREdistribution['group B']={}

for gene in GREgenes:

    if reverseAnnotationMap[gene] in geneSets['group A']:
        groupLabel='group A'
    elif reverseAnnotationMap[gene] in geneSets['group B']:
        groupLabel='group B'
    else:
        print('error while associating gene to a group')
        sys.exit()

    if '9.and.28' not in GREdistribution[groupLabel]:
        GREdistribution[groupLabel]['9.and.28']=0
    if '7' not in GREdistribution[groupLabel]:
        GREdistribution[groupLabel]['7']=0
    if 'others' not in GREdistribution[groupLabel]:
        GREdistribution[groupLabel]['others']=0

    for dimension in dimensions:
        if dimension in geneGREs[gene]:
            value=geneGREs[gene][dimension]
        else:
            value=0.
            
        # adding to GREdistribution
        if dimension == 9 or dimension == 28:
            GREdistribution[groupLabel]['9.and.28']=GREdistribution[groupLabel]['9.and.28']+value
        elif dimension == 7:
            GREdistribution[groupLabel]['7']=GREdistribution[groupLabel]['7']+value
        else:
            GREdistribution[groupLabel]['others']=GREdistribution[groupLabel]['others']+value
            
print(GREdistribution)                                




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
