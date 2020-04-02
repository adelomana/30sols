import sys,numpy,scipy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.cluster,sklearn.datasets

# 0. define user variables
print('initializing...')
dataFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.txt'
outputFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.clean.csv'
ribosomalGenesFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
ribosomal_annotation_file = '/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

# 1. read data
print('reading data...')

# 1.1. read ribosomal gene names
RP_gene_names = []
with open(ribosomal_annotation_file, 'r') as f:
    next(f)
    for line in f:
        vector = line.split('\t')
        vng = vector[1]
        RP_gene_names.append(vng)

# 1.2. reading corem motif counts
data={} # data[coremName]=[(GRE,frequency),(GRE,frequency),...]

allCoremNames=[]
allGREnames=[]
all_non_RPs = []

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

        # filter out corems that have < 33% RPs and > 2 non-RPs
        RP_members = []; non_RP_members = []
        for element in gene_memberships:
            if element in RP_gene_names:
                RP_members.append(element)
            else:
                non_RP_members.append(element)
            
        co = len(gene_memberships)
        rb = len(RP_members)
        nrb = len(non_RP_members)
        purity = rb / co

        for element in non_RP_members:
            if element not in all_non_RPs:
                all_non_RPs.append(element)

        if rb > 2 and purity > 1/3:
            print('purity', purity)
            print('ribo : non ribo', rb, nrb)
            print()
        
            if coremName not in data.keys():
                data[coremName]={}

            if motifName not in data[coremName].keys():
                data[coremName][motifName]=frequency
        
            if coremName not in allCoremNames:
                allCoremNames.append(coremName)
            if motifName not in allGREnames:
                allGREnames.append(motifName)

# prepare file for enrichment. Strange format to make it compatible
print('non RP genes', len(all_non_RPs))
file_name = '/Users/alomana/github/30sol/F1.interplay/panel.a/results/results.others.txt'
with open(file_name, 'w') as f:
    for element in all_non_RPs:
        f.write('{}\t{}\t{}\n'.format('other', element, 'tempo'))
print('done')

# 2. arranging data
print('preparing data...')
coremInt=[int(element) for element in allCoremNames]
coremInt.sort()
sortedCoremNames=[str(element) for element in coremInt]

GREnamesInt=[int(element.split('_')[1]) for element in allGREnames]
GREnamesInt.sort()
sortedGREnames=['MOTC_{}'.format(element) for element in GREnamesInt]

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
