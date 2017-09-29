import sys,numpy,scipy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.cluster,sklearn.datasets


# 0. define user variables
print('initializing...')
dataFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.txt'
outputFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.clean.csv'

# 1. reading data
print('reading data...')
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
        
        if coremName not in data.keys():
            data[coremName]={}

        if motifName not in data[coremName].keys():
            data[coremName][motifName]=frequency
        
        if coremName not in allCoremNames:
            allCoremNames.append(coremName)
        if motifName not in allGREnames:
            allGREnames.append(motifName)

# 2. arranging data
print('preparing data...')
coremInt=[int(element) for element in allCoremNames]
coremInt.sort()
sortedCoremNames=[str(element) for element in coremInt]

GREnamesInt=[int(element.split('_')[1]) for element in allGREnames]
GREnamesInt.sort()
sortedGREnames=['MOTC_{}'.format(element) for element in GREnamesInt]


# writing the output file
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
