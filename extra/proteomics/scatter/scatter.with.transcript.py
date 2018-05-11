###
### This script creates a scatter plot of log2 FC of transcript and pt abundance for ribo-pts
###

import os,sys,numpy,seaborn,pandas
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def dataReader():

    '''
    This function reads data available and outputs the defined dictionary.
    '''

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC

    conditions=[]
    geneNames=[]
    timepoints=[]
    replicates=[]
    
    allFiles=os.listdir(dataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=dataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in conditions:
            conditions.append(condition)
        if replicate not in replicates:
            replicates.append(replicate)

        if condition not in data.keys():
            data[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}

        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')

                geneName=vector[0]
                if geneName not in geneNames:
                    geneNames.append(geneName)

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c

    # sort
    geneNames.sort()
    conditions.sort()
    replicates.sort()
    
    return data,geneNames,conditions,replicates,timepoints

###
### MAIN
###

# 0. user defined variables
dataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

# 1. read data
print('reading data...')
riboPtNames=riboPtNamesReader()
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteinDataReader()

# 2. analysis
