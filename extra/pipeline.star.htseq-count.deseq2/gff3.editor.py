import sys

inputFileName='/Volumes/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.gff3'
outputFileName='/Volumes/omics4tb/alomana/projects/TLR/data/genome/hsa.ASM680v1.edited.gff3'

# 1. learning the equivalences of gene0 AE004437.1 to final gene name
mapOfids={}
allGeneNames=[]
with open(inputFileName,'r') as f:
    for line in f:
        vector=line.split('\t')
        if vector[0][0] != '#':
            if vector[2] == 'gene':
                contig=vector[0]
                info=vector[-1]
                rankName=info.split('ID=')[1].split(';')[0]
                
                if contig == 'AE004437.1':
                    geneName=info.split('locus_tag=')[1].split(';')[0]
                    geneName=geneName.replace('\n','')
                elif contig == 'AF016485.1':
                    geneName=info.split('Name=')[1].split(';')[0]
                    if geneName in ['sojB','tbpC','repI']: # three genes are defined at both plasmids
                        geneName=geneName+'.AF016485'
                elif contig == 'AE004438.1':
                    geneName=info.split('Name=')[1].split(';')[0]
                else:
                    print('contig not found')
                    sys.exit()

                label=contig+'.'+rankName
                mapOfids[contig+'.'+rankName]=geneName
                allGeneNames.append(geneName)

# 2. checking for non unique names
uniques=list(set(allGeneNames))
if len(allGeneNames) != len(uniques):
    print('repeated annotation')
    sys.exit()

# 2. substituting the values
g=open(outputFileName,'w')

with open(inputFileName,'r') as f:
    for line in f:
        vector=line.split('\t')
        if vector[0][0] != '#':
            contig=vector[0]
            info=vector[-1]


            if vector[2] == 'gene' or vector[2] == 'CDS':
                info=vector[-1]
                if '=gene' in info:
                    
                    gene2mod='gene'+info.split('=gene')[1].split(';')[0]
                    dictLabel=contig+'.'+gene2mod
                    newLine=line.replace(gene2mod,mapOfids[dictLabel])
                    
            else:
                newLine=line
        else:
            newLine=line
        g.write(newLine)
g.close()
