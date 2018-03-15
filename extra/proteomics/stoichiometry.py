###
### This script computes the stoichiometry of ribosomal proteins over time
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

def figureGrapher():

    '''
    This function creates a figure on the stoichiometry analysis.
    '''

    print('\nanalysis of condition {}...'.format(condition))
    
    # f.1. initialize variables
    foldChangeInfo={} # dictionary with all information about fold-change: foldChangeInfo[timePointLabel][riboPtName]=value (fold-change)
    stoichInfo={} # dictionary with all information about stoichiometry: stoichInfo[timePointLabel][riboPtName]=value (log2 stoichiometry value)

    timeStampsViolin=[]; stoichValuesViolin=[]
    timeStampsSwarm=[]; stoichValuesSwarm=[]
    timeStampsViolin.append(timePointLabels[0]); stoichValuesViolin.append(0) # incorporate a single point for violin
    timeStampsSwarm.append(timePointLabels[0]); stoichValuesSwarm.append(0) # incorporate a single point for violin

    significantCases=[]
    significantPositions=[]
    noInfoSet=[]

    # f.2. compute stoichiometries
    for i in range(len(timepoints)):

        timeLabel=timePointLabels[i+1]
        
        if timeLabel not in foldChangeInfo:
            foldChangeInfo[timeLabel]={}

        if timeLabel not in stoichInfo:
            stoichInfo[timeLabel]={}
            
        for ribopt in riboPtNames:
            try:
                x=numpy.median([data[condition][replicate][timepoints[i]][ribopt] for replicate in replicates])
                foldChangeInfo[timeLabel][ribopt]=2**x
            except:
                if ribopt not in noInfoSet:
                    noInfoSet.append(ribopt)
                    print('no data for {} {}'.format(timepoints[i],ribopt))

        # f.2.1. compute the stoichiometry per time point
        localNames=list(foldChangeInfo[timeLabel].keys())
        allFractions=[foldChangeInfo[timeLabel][localName] for localName in localNames]
        theSum=sum(allFractions)
        stoich=(numpy.array(allFractions)/theSum)*len(allFractions)
        log2Stoich=numpy.log2(stoich)
        for j in range(len(localNames)):
            stoichInfo[timeLabel][localNames[j]]=log2Stoich[j]

        # f.2.2. finding the limits of 95% of the distribution
        zp=1.959963984540 # taken from https://en.wikipedia.org/wiki/Normal_distribution
        mean=numpy.mean(log2Stoich)
        standardDeviation=numpy.std(log2Stoich)
        low=mean-standardDeviation; high=mean+standardDeviation
        
        print('low, high',low,high)

        # f.2.3. fill up variables for plotting considering significances
        for i in range(len(localNames)):
            v=log2Stoich[i]
            timeStampsViolin.append(timeLabel); stoichValuesViolin.append(v)
            if v > high or v < low:
                significantPositions.append([timePointLabels.index(timeLabel),v])
                if localNames[i] not in significantCases:
                    significantCases.append(localNames[i])
            else:
                timeStampsSwarm.append(timeLabel); stoichValuesSwarm.append(v)

    # f.5. create a dataframe for plotting with seaborn
    stoichiometryViolin=list(zip(timeStampsViolin,stoichValuesViolin))
    dfViolin=pandas.DataFrame(data=stoichiometryViolin,columns=['Time points','Stoichiometry'])

    stoichiometrySwarm=list(zip(timeStampsSwarm,stoichValuesSwarm))
    dfSwarm=pandas.DataFrame(data=stoichiometrySwarm,columns=['Time points','Stoichiometry'])
    
    # f.6. plot violin and swarm plots with seaborn
    colorPalette={}
    for i in range(len(orderedColors)):
        colorPalette[timePointLabels[i]]=orderedColors[i]

    ax=seaborn.violinplot(x='Time points',y='Stoichiometry',data=dfViolin,inner=None,linewidth=0,color='0.5')
    matplotlib.pyplot.setp(ax.collections, alpha=.5)

    ax=seaborn.swarmplot(x='Time points',y='Stoichiometry',data=dfSwarm,color='white',size=theDotSize,zorder=2)

    # f.7. plot special points
    matplotlib.pyplot.plot(0,0,'o',color=orderedColors[0],ms=theDotSize,mew=0,zorder=10)

    # f.8. plot significant cases
    for case in significantCases:
        x=[0]
        y=[0]
        for i in range(len(timepoints)):
            timeLabel=timePointLabels[i+1]
            x.append(i+1)
            y.append(stoichInfo[timeLabel][case])
            matplotlib.pyplot.plot(i+1,stoichInfo[timeLabel][case],'o',color='magenta',ms=theDotSize,mew=0)
        matplotlib.pyplot.plot(x,y,':',color='black',lw=1,alpha=0.3,zorder=0)

    for position in significantPositions:
        matplotlib.pyplot.plot(position[0],position[1],'sk',ms=theDotSize*2,mew=0)

    # f.7. final figure closing
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Stoichiometry (log$_2$ ribo-pt)')
    matplotlib.pyplot.title(condition)

    #matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,ncol=1,fontsize=14)

    figureName='figure.{}.pdf'.format(condition)
    matplotlib.pyplot.tight_layout()
    #matplotlib.pyplot.axes().set_aspect('equal')
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

orderedColors=['red','orange','green','blue']
timePointLabels=['TP1','TP2','TP3','TP4']

theDotSize=3

conditions=['rbf','lysate']

# 1. read data
print('reading data...')
data,geneNames,conditions,replicates,timepoints=dataReader()
riboPtNames=riboPtNamesReader()

# 2. analyse data
print('analyzing data...')
for condition in conditions:
    figureGrapher()
