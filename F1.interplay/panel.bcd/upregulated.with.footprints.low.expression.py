#transcripts upregulated but have higher footprints than expected, are they at low abundance?


#define on the fc figure, the deviation from PIs.
#Then for those who are upreguated, those that have greater changes in



###

#either, run PI with all and then find if tc plus footprints have lower exp than tc plus all

#or
#get the names of tc plus
#get the names of high or low footprints
#are low and high different levels of exp???

###
### This script tests the question: what is the transcript abundance difference in TP1 between transcripts with more or less footprints than expected among the upregulated transcripts?
###

import sys

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    geneNames=[]; timepoints=[]; replicates=[]
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

            if replicate not in replicates:
                replicates.append(replicate)
            if timepoint not in timepoints:
                timepoints.append(timepoint)

            if fraction not in data.keys():
                data[fraction]={}
            if replicate not in data[fraction].keys():
                data[fraction][replicate]={}
            if timepoint not in data[fraction][replicate].keys():
                data[fraction][replicate][timepoint]={}
            
        for line in f:
            vector=line.split('\t')[:-1]
            values=[float(element) for element in vector[1:]]
            geneName=vector[0].replace('_','')
            if geneName not in geneNames:
                geneNames.append(geneName)
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data,geneNames,timepoints,replicates

# 0. user defined variables

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'
deviationFile='diagonal.deviated.names.txt'
upregulatedTranscriptsFile='../panel.a/results/results.black.plus.txt'

# 1. read data

# 1.1. retrieve up-regulated transcripts
upregulatedTranscripts=[]
with open(upregulatedTranscriptsFile,'r') as f:
    for line in f:
        v=line.split('\t')
        name=v[2].replace('gene-','')
        upregulatedTranscripts.append(name)

# 1.2. read up, neutral or down for footprints
diagonalDeviations={} # diagonalDeviations[TP1|TP2|TP3|TP4][up|neutral|down]=[name1,name2,name3,...]
with open(deviationFile,'r') as f:
    for line in f:
        v=line.split()
        timepoint=v[0]
        regulation=v[1]
        name=v[2]

        if timepoint not in diagonalDeviations:
            diagonalDeviations[timepoint]={}
        if regulation not in diagonalDeviations[timepoint]:
            diagonalDeviations[timepoint][regulation]=[]
        if name not in diagonalDeviations[timepoint][regulation]:
            diagonalDeviations[timepoint][regulation].append(name)
                
# 1.3. reading mRNA data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 2. process data

#select genes that have specified footprints
timepoints=list(diagonalDeviations.keys())
fragileGroups={}
fragileGroups['up']=diagonalDeviations[timepoints[0]]['up']
fragileGroups['neutral']=diagonalDeviations[timepoints[0]]['neutral']
fragileGroups['down']=diagonalDeviations[timepoints[0]]['down']

robustGroups={}
robustGroups['up']=[]; robustGroups['neutral']=[]; robustGroups['down']=[]
up=[]; neutral=[]; down=[]
for element in diagonalDeviations[timepoints[0]]['up']:
    successes=0
    for timepoint in timepoints:
        if element in diagonalDeviations[timepoint]['up']:
            successes=successes+1
    if successes == 4:
        robustGroups['up'].append(element)

for element in diagonalDeviations[timepoints[0]]['neutral']:
    successes=0
    for timepoint in timepoints:
        if element in diagonalDeviations[timepoint]['neutral']:
            successes=successes+1
    if successes == 4:
        robustGroups['neutral'].append(element)

for element in diagonalDeviations[timepoints[0]]['down']:
    successes=0
    for timepoint in timepoints:
        if element in diagonalDeviations[timepoint]['down']:
            successes=successes+1
    if successes == 4:
        robustGroups['down'].append(element)

# intersect up and down with upregulated genes
a=list(set(upregulatedTranscripts) & set(fragileGroups['up']))
b=list(set(upregulatedTranscripts) & set(fragileGroups['down']))

c=list(set(upregulatedTranscripts) & set(robustGroups['up']))
d=list(set(upregulatedTranscripts) & set(robustGroups['down']))

print(len(a),len(b),len(c),len(d))

print([rnaExpression['trna']['br1']['tp.1'][element] for element in b])

print(len(upregulatedTranscripts),len(fragileGroups['up']))
#make a figure with three histograms and run a hypothesis test


# build a figure that shows expression of upregulated genes into three groups: consistent footprints up, neutral or down
# check first in TP1, then cosistently across all TPs
#print(DETs['tp.4'].keys())

