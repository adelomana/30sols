###
### This script computes the stoichiometry of ribosomal proteins over time
###

import os,sys,numpy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.decomposition,sklearn.manifold

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

def figureGrapher():

    '''
    This function creates a figure on the stoichiometry analysis.
    '''

    print(condition)
    
    # f.1. plotting the first point
    matplotlib.pyplot.plot([1],[1],'o',color=colors[0],label=timePointLabels[0],alpha=0.3,mew=0,ms=theDotSize)

    # f.2. processing 
    noInfoSet=[]
    for i in range(len(timepoints)):
        allFractions=[]
        for ribopt in riboPtNames:
            try:
                x=numpy.median([data[condition][replicate][timepoints[i]][ribopt] for replicate in replicates])
                allFractions.append(2**x)
            except:
                if ribopt not in noInfoSet:
                    noInfoSet.append(ribopt)
                    print('no data for {} {}'.format(timepoints[i],ribopt))

        # f.2. computing the stoichiometry per time point
        theSum=sum(allFractions)
        stoich=(numpy.array(allFractions)/theSum)*len(allFractions)

        print(allFractions)
        print(theSum)
        print(stoich)
        #sys.exit()

        matplotlib.pyplot.plot(numpy.repeat(i+2,len(stoich)),stoich,'o',color=colors[i+1],alpha=0.3,mew=0,ms=theDotSize,label=timePointLabels[i+1])

    # do consistently the top and bottome 5 percentile
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Ribosomal stoichiometry')
    matplotlib.pyplot.title(condition)

    matplotlib.pyplot.xticks([1,2,3,4])
    #matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,ncol=1,fontsize=14)

    figureName='figure.{}.pdf'.format(condition)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()    

    return None

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
### MAIN
###

# 0. user defined variables
dataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

colors=['red','orange','green','blue']
timePointLabels=['TP1','TP2','TP3','TP4']
theDotSize=8

conditions=['rbf','lysate']

# 1. read data
print('reading data...')
data,geneNames,conditions,replicates,timepoints=dataReader()
riboPtNames=riboPtNamesReader()
print(timepoints)

# 2. analyse data
print('analyzing data...')
for condition in conditions:
    figureGrapher()
