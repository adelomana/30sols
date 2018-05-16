###
### This script investigates if protein changes can be explained by ribosomal footprints, i.e., translational regulation.
###

import numpy,sys,os

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    This function creates a dictionary for expression values as
    expression[trna/rbf][ribo-pt gene name][timepoint][replicate]=value
    '''

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    with open(expressionDataFile,'r') as f:

        firstLine=f.readline()
        header=firstLine.split(',')
        sampleNames=header[1:]
        sampleNames[-1]=sampleNames[-1].replace('\n','')

        for line in f:
            vector=line.split(',')

            # geneName
            geneName=vector[0]
            if geneName not in geneNames:
                geneNames.append(geneName)

            for i in range(len(sampleNames)):

                # sampleType
                sampleType=sampleNames[i].split('.')[0]
                if sampleType not in sampleTypes:
                    sampleTypes.append(sampleType)

                # timepoint
                timepoint='tp.{}'.format(int(sampleNames[i].split('.')[-1]))
                if timepoint not in timepoints:
                    timepoints.append(timepoint)

                # replicate
                replicate='rep.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
                if replicate not in replicates:
                    replicates.append(replicate)

                # value
                value=float(vector[i+1])

                # make sure keys exist
                if sampleType not in expression.keys():
                    expression[sampleType]={}
                if geneName not in expression[sampleType].keys():
                    expression[sampleType][geneName]={}
                if timepoint not in expression[sampleType][geneName].keys():
                    expression[sampleType][geneName][timepoint]={}

                expression[sampleType][geneName][timepoint][replicate]=value

    # sort variables
    sampleTypes.sort()
    geneNames.sort()
    timepoints.sort()
    replicates.sort()

    return expression,sampleTypes,geneNames,timepoints,replicates

def proteomicsDataReader():

    '''
    This function reads data available and outputs the defined dictionary.
    '''

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC

    conditions=[]
    geneNames=[]
    timepoints=[]
    replicates=[]
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

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
                consider here converting old to new names
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

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])
            
    return riboPtNames

###
### M A I N
###

# 0. user defined variables.
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'

# 1. read data.
print('reading data...')
riboPtNames=riboPtNamesReader()
print('\t reading transcriptome data...')
expression,expressionSampleTypes,expressionGeneNames,expressionTimepoints,expressionReplicates=expressionReader()

#print(expression)
print(expressionSampleTypes)
print(expressionGeneNames[:10])
print(expressionTimepoints)
print(expressionReplicates)

print('\t reading proteomics data...')
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteomicsDataReader()
print(proteinNames[:10])
print(proteinConditions)
print(proteinReplicates)
print(proteinTimepoints)

# 2. filter data.
# Filter data. Remove pt/mRNA/RF whose standard error of the mean is larger than 30%.
print('filtering data...')

# 3. analysis
print('running analysis...')

# 3.1. scatter plot of FC_pt vs FC_mRNA. This figure would reveal how well mRNA explains pt changes.
print('\t building pt vs mRNA...')

# x.2. scatter plot of FC_pt vs FC_RF. This figure would reveal how well RF explains pt changes.
print('\t building pt vs RF...')

# x.3 scatter plot of FC_pt vs FC_mRNA + FC_RF. This figure would reveal how well an equal linear model would explain pt changes.
print('\t building pt vs equal model...')

# x.4. scatter plot of model: FC_pt = w_1 FC_mRNA + w_2 (FC_RF - FC_RF_0).
# For this section I need to obtain FC_RF_0 from linear regression on FC_mRNA vs FC_RF, from all genes.
# w_1 is defined within the range [0,1].
# w_2 is defined within the range [-1,1].
print('\t building weighted model...')

# x.5. visualization of w_1 vs w_2 to distinguish TCR from TLR.
# I should have one for each time point - t0
print('\t visualizing regulatory weights...')

# x. a sanity check for the script is the capture of RF causing pausing on rps10-like.
