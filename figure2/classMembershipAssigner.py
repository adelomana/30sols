###
### this script assigns most voted class to corems
###

import sys

# 0. user defined variables
gene2ClassFile='/Volumes/omics4tb/alomana/projects/TLR/data/annotation/si.table.1.si.ribosomal.protein.index.information.csv'
geneDictionaryFile='/Volumes/omics4tb/alomana/projects/TLR/data/annotation/si.table.1.gene.dictionary.csv'
coremsFile='/Volumes/omics4tb/alomana/projects/TLR/data/annotation/si.table.1.corem.index.information.txt'

outputFile='/Volumes/omics4tb/alomana/projects/TLR/data/annotation/coremClasses.csv'

acceptedClasses=['Blue','Red','Green','Magenta']

# 1. building a dictionary of names
synonyms={}
with open(geneDictionaryFile,'r') as f:
    for line in f:
        vector=line.split(',')
        geneName=vector[1]
        riboName=vector[2].replace('\n','')
        synonyms[geneName]=riboName

# 2. associate genes to class
gene2Class={}
with open(gene2ClassFile,'r') as f:
    for line in f:
        vector=line.split(',')
        name=vector[1]
        theClass=vector[2]
        if theClass in acceptedClasses:
            gene2Class[synonyms[name]]=theClass

# 3. associate corems to class
g=open(outputFile,'w')
g.write('coremID,ribosomal.class\n')
coremClasses={}
with open(coremsFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split('\t')
        coremID=int(vector[1])
        pre=vector[4]
        pro=pre.split(',')
        pru=[element.replace('"','') for element in pro]
        proteins=[element.replace(' ','') for element in pru]

        foundClasses=[gene2Class[element] for element in proteins]
        democracy=max(set(foundClasses), key=foundClasses.count)
        
        g.write('corem.{},{}\n'.format(coremID,democracy.lower()))

        print(democracy)
        
g.close()
