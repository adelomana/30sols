"""
this script finds which other genes are present in the ribosomal corems that are not ribosomal genes
"""

import sys

# 0. user defined data
corems2GenesFile='/Volumes/omics4tb/alomana/projects/TLR/data/corem/hsa_c2g.txt'
ribosomalGenesFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
ribosomalCoremsFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/classMembership.txt'

# 1. reading files

"""
select ribosomal corems
find all genes
substract  ribosomal genes
"""

# 1.1. reading ribosomal gene names
ribosomalGenes=[]
with open(ribosomalGenesFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split()
        ribosomalGenes.append(vector[0])

# 1.2. reading ribosomal corems
ribosomalCorems=[]
with open(ribosomalCoremsFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split()
        ribosomalCorems.append(vector[1])

# 1.3. reading all genes of corems
allGenes=[]
with open(corems2GenesFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split()
        corem=vector[0]
        if corem in ribosomalCorems:
            some=vector[1].split(';')
            print(some)
            for element in some:
                if element not in allGenes:
                    allGenes.append(element)

# 2. analysis

# 2.1. check all ribosomal genes are recapitulated
intersect=list(set(ribosomalGenes) & set(allGenes))
print(len(intersect),len(ribosomalGenes))

# 2.2. find the set of genes that are not ribosomal genes
corregulated=[]
for element in allGenes:
    if element not in ribosomalGenes:
        corregulated.append(element)

print(corregulated,len(corregulated))
