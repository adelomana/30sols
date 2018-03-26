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
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
operonPredictionsFile='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/operonPredictions.txt'
operonPredictionsDir='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/'

# 1. read data

# 1.1. ribosomal genes
riboPtNames,synonyms=riboPtNamesReader()

# 1.2. operon predictions
operonPredictions=operonPredictionsReader()

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

# 4. saving files

# 4.1. ribo-pt gene operons

