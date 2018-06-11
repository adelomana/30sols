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

def figureGrapher(colorAssociation):

    '''
    This function creates a figure on the stoichiometry analysis.
    '''

    print('\nanalysis of condition {}...'.format(condition))
    
    # f.1. initialize variables
    noInfoSet=[]
    
    foldChangeInfo={} # dictionary with all information about fold-change: foldChangeInfo[timePointLabel][riboPtName]=value (fold-change)
    stoichInfo={} # dictionary with all information about stoichiometry: stoichInfo[timePointLabel][riboPtName]=value (log2 stoichiometry value)

    timeStampsViolin=[]; stoichValuesViolin=[]
    timeStampsSwarm=[]; stoichValuesSwarm=[]
    timeStampsViolin.append(timePointLabels[0]); stoichValuesViolin.append(0) # incorporate a single point for violin
    timeStampsSwarm.append(timePointLabels[0]); stoichValuesSwarm.append(0) # incorporate a single point for violin

    significantNames=[]
    significantPositions={}

    thresholds={}

    # f.2. compute stoichiometries
    for i in range(len(timepoints)):

        timeLabel=timePointLabels[i+1]
        
        if timeLabel not in foldChangeInfo:
            foldChangeInfo[timeLabel]={}

        if timeLabel not in stoichInfo:
            stoichInfo[timeLabel]={}
            
        for ribopt in riboPtNames:
            values=[]
            for replicate in replicates:
                value=None
                try:
                    value=data[condition][replicate][timepoints[i]][ribopt]
                except:
                    pass
                if value != None:
                    values.append(value)
            if len(values) >= 3 :
                average=numpy.median(values)
                sem=numpy.std(values)/numpy.sqrt(len(values))
                rsem=sem/numpy.mean(values)
                if rsem < 0.3: 
                    foldChangeInfo[timeLabel][ribopt]=2**average
                else:
                    print('\t\t\t loosing {} for low precision: {} {}'.format(ribopt,values,rsem))
            else:
                 print('\t\t loosing {} for not enough replicates: {}'.format(ribopt,values))                

        # f.2.1. compute the stoichiometry per time point
        localNames=list(foldChangeInfo[timeLabel].keys())
        print('{} {} n = {}'.format(timeLabel,condition,len(localNames)))
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
        
        margin=0.2
        a=i+1-margin
        b=i+1+margin
        matplotlib.pyplot.plot([a,b],[high,high],'-',color='white')
        matplotlib.pyplot.plot([a,b],[low,low],'-',color='white')

        # f.2.3. fill up variables for plotting considering significances
        for i in range(len(localNames)):
            v=log2Stoich[i]
            timeStampsViolin.append(timeLabel); stoichValuesViolin.append(v)
            if v > high or v < low:

                if localNames[i] not in significantPositions:
                    significantPositions[localNames[i]]=[[],[]]
                significantPositions[localNames[i]][0].append(timePointLabels.index(timeLabel))
                significantPositions[localNames[i]][1].append(v)
                
                if localNames[i] not in significantNames:
                    significantNames.append(localNames[i])
            else:
                timeStampsSwarm.append(timeLabel); stoichValuesSwarm.append(v)

    # f.3. create a dataframe for plotting with seaborn
    stoichiometryViolin=list(zip(timeStampsViolin,stoichValuesViolin))
    dfViolin=pandas.DataFrame(data=stoichiometryViolin,columns=['Time points','Stoichiometry'])

    stoichiometrySwarm=list(zip(timeStampsSwarm,stoichValuesSwarm))
    dfSwarm=pandas.DataFrame(data=stoichiometrySwarm,columns=['Time points','Stoichiometry'])
    
    # f.4. plot violin and swarm plots with seaborn
    ax=seaborn.violinplot(x='Time points',y='Stoichiometry',data=dfViolin,inner=None,linewidth=0,color='0.5')
    matplotlib.pyplot.setp(ax.collections, alpha=.5)
    ax=seaborn.swarmplot(x='Time points',y='Stoichiometry',data=dfSwarm,color='white',size=theDotSize,zorder=1)

    # f.5. plot special point
    matplotlib.pyplot.plot(0,0,'s',color='black',ms=theDotSize,mew=0,zorder=10)

    # f.6. plot significant trajectories
    for name in sorted(significantPositions):
        if len(significantPositions[name][0]) > 1:
            x=[0]; y=[0]
            for i in range(len(timepoints)):
                timeLabel=timePointLabels[i+1]
                x.append(i+1)
                y.append(stoichInfo[timeLabel][name])
            if name not in colorAssociation:
                colorAssociation[name]=matplotlib.cm.tab10(len(colorAssociation))
            matplotlib.pyplot.plot(x,y,':',color=colorAssociation[name],lw=3,zorder=0,label=nameAliases[name])

    # f.7. plot significant points
    for name in sorted(significantPositions):
        x=significantPositions[name][0]
        y=significantPositions[name][1]
        if name in colorAssociation:
            matplotlib.pyplot.plot(x,y,'o',color=colorAssociation[name],ms=theDotSize*2.5,mew=0)
        else:
            matplotlib.pyplot.plot(x,y,'o',color='black',ms=theDotSize,mew=0)
            #matplotlib.pyplot.text(x[0],y[0],nameAliases[name])
            
    # f.8. final figure closing
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Stoichiometry (log$_2$ ribo-pt)')
    matplotlib.pyplot.title(condition)

    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=3,ncol=2,fontsize=14)

    figureName='figure.{}.pdf'.format(condition)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()    

    return colorAssociation

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[1])
            
    return riboPtNames

###
### MAIN
###

# 0. user defined variables
dataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

timePointLabels=['TP1','TP2','TP3','TP4']
theDotSize=3
conditions=['rbf','lysate']

nameAliases={}
nameAliases['VNG0551G']='L44E'
nameAliases['VNG1158G']='S28E'
nameAliases['VNG2469G']='L39E'
nameAliases['VNG0433C']='S10-like'
nameAliases['VNG1132G']='S13'
nameAliases['VNG1159G']='L24E'

nameAliases['VNG1705G']='L5'
nameAliases['VNG1695G']='L22'
nameAliases['VNG1703G']='S4E'
nameAliases['VNG1697G']='S3'
nameAliases['VNG2048G']='S24E'
nameAliases['VNG1134G']='S11'
nameAliases['VNG1157G']='L7AE'
nameAliases['VNG2047G']='S27AE'

# 1. read data
print('reading data...')
data,geneNames,conditions,replicates,timepoints=dataReader()
riboPtNames=riboPtNamesReader()

# 2. analyse data
print('analyzing data...')
colorAssociation={}
for condition in conditions:
    colorAssociation=figureGrapher(colorAssociation)
