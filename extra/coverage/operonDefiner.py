###
### This script formats ribo-pt operons from MicrobesOnline, amenable to coverage analysis.
###

import sys

def operonPredictionsReader():

    '''
    This function reads operon predictions from Microbes Online: http://www.microbesonline.org/operons/gnc64091.html
    '''

    operonPredictions={}
    operon=[]
    count=0
    
    with open(operonPredictionsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')

            gene1=vector[2]
            gene2=vector[3]
            prediction=vector[6]

            if prediction == 'TRUE':

                if operon == []:
                    operon=[gene1,gene2]
                else:
                    if gene1 in operon:
                        operon.append(gene2)
            else:
                if operon != []:
                    count=count+1
                    operonPredictions['OP{}'.format(str(count).zfill(4))]=operon
                    operon=[]

    return operonPredictions

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[1])

    riboPtNames.sort()
            
    return riboPtNames

def synonymsReader():

    '''
    This function reads the GFF3 file and returns a dictionary with synonyms between old and new locus names.
    '''

    synonyms={}
    with open(annotationFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1].replace('\n','')
                if 'old_locus_tag=' in info:
                    old=info.split('old_locus_tag=')[1].split(';')[0]
                    new=info.split('ID=')[1].split(';')[0]

                    if '%' in old:
                        olds=old.split('%2C')
                        for element in olds:
                            synonyms[element]=new
                    else:
                        synonyms[old]=new
                
    return synonyms

###
### MAIN
###

# 0. user defined variables
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
operonPredictionsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/operonPredictions.txt'
operonPredictionsDir='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

# 1. read data

# 1.1. ribosomal genes
riboPtNames=riboPtNamesReader()

# 1.2. operon predictions
operonPredictions=operonPredictionsReader()

# 1.3. synonyms definer
synonyms=synonymsReader()

# 2. save ribo-pt operons containing ribo-pt genes
riboOperons=[]
for operon in sorted(operonPredictions):
    localGenes=operonPredictions[operon]
    for name in localGenes:
        if name in riboPtNames:
            if operon not in riboOperons:
                riboOperons.append(operon)

# 3. check that all ribo-pt genes are in the ribo-operons
print('Ribo-pt genes found: {}'.format(len(riboPtNames)))

found=[]
for operon in riboOperons:
    for element in operonPredictions[operon]:
        if element in riboPtNames:
            if element not in found:
                found.append(element)
                
print('Rib-pt genes found in operons: {}'.format(len(found)))

# 3.1. define the non-operon ribo-pt genes (NORPGs)
NORPGs=[]
for name in riboPtNames:
    if name not in found:
        NORPGs.append(name)
NORPGs.sort()

# 4. saving files

# 4.1. ribo-pt gene operons
fileName=operonPredictionsDir+'riboPtOperons.txt'
with open(fileName,'w') as f:
    f.write('# operonID\tgeneNames\n')
    for operon in riboOperons:
        f.write('{}'.format(operon))
        for name in operonPredictions[operon]:
            f.write('\t{}'.format(synonyms[name]))
        f.write('\n')

# 4.2. non-operon ribo-pt genes
fileName=operonPredictionsDir+'NORPGs.txt'
with open(fileName,'w') as f:
    f.write('# This file contains non-operon ribo-pt genes (NORPGs)\n')
    for name in NORPGs:
        f.write('{}\n'.format(synonyms[name]))
