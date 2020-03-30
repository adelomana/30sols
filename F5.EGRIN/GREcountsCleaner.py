import sys,numpy,scipy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.cluster,sklearn.datasets


# 0. define user variables
print('initializing...')
dataFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.txt'
outputFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.clean.csv'
ribosomalGenesFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

# 1. read data
print('reading data...')

# 1.1. read ribosomal gene names
ribosomalGenes=[]
with open(ribosomalGenesFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split()
        ribosomalGenes.append(vector[0])
# convert RP gene IDs to new annotation


# 1.2. reading corem motif counts
data={} # data[coremName]=[(GRE,frequency),(GRE,frequency),...]

allCoremNames=[]
allGREnames=[]

with open(dataFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split('\t')

        coremName=vector[0].replace('"','')
        motifName=vector[1].replace('"','')
        frequencyString=vector[2].replace('"','')
        frequency=int(frequencyString)
        gene_memberships = vector[3].split(':')
        gene_memberships = [element.replace('"', '') for element in gene_memberships]
        gene_memberships = [element.replace('\n', '') for element in gene_memberships]

        is_RP_corem = False # filter out corems that have < 80% RPs and > 5 non-RPs
        labels = [element in ribosomalGenes for element in gene_memberships]
        print(gene_memberships[:10])
        print(ribosomalGenes[:10])
        print(labels)
        sys.exit()
        

        print(coremName, motifName, frequency)
        print(gene_membership)
        print()


        
        if coremName not in data.keys():
            data[coremName]={}

        if motifName not in data[coremName].keys():
            data[coremName][motifName]=frequency
        
        if coremName not in allCoremNames:
            allCoremNames.append(coremName)
        if motifName not in allGREnames:
            allGREnames.append(motifName)

sys.exit()

# 2. arranging data
print('preparing data...')
coremInt=[int(element) for element in allCoremNames]
coremInt.sort()
sortedCoremNames=[str(element) for element in coremInt]

GREnamesInt=[int(element.split('_')[1]) for element in allGREnames]
GREnamesInt.sort()
sortedGREnames=['MOTC_{}'.format(element) for element in GREnamesInt]

# 3. 
for corem in sortedCoremNames:
    print(corem)
sys.exit()


# 4. writing the output file
g=open(outputFile,'w')

g.write('coremName')
for gre in sortedGREnames:
    g.write(',{}'.format(gre))
g.write('\n')

for corem in sortedCoremNames:
    g.write('corem.{}'.format(corem))    
    for gre in sortedGREnames:
        try:
            value=data[corem][gre]
        except:
            value=0
        g.write(',{}'.format(value))
    g.write('\n')

g.close()
