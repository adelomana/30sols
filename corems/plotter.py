import os,sys,numpy,copy
import matplotlib,matplotlib.pyplot,matplotlib.cm,matplotlib.patches

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def colorDefiner(sortedConditions):

    # mapping annotation to values
    values=[annotationValues.index(sampleMetadata[element]) for element in sortedConditions]

    # mapping values to colors
    customMap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min(values),vmax=max(values)),cmap='tab20')
    theColors=[customMap.to_rgba(value) for value in values]

    # building a dictionary between colors and annotation for legend
    legendColors={}
    for i in range(len(sortedConditions)):
        if sampleMetadata[sortedConditions[i]] not in legendColors.keys():
            legendColors[sampleMetadata[sortedConditions[i]]]=theColors[i]

    return theColors,legendColors

def coremReader(inputFileName):

    expression=[]

    with open(inputFileName,'r') as f:
        
        header=f.readline()
        vector=header.split('\t')
        conditions=[element.replace('"','') for element in vector]
        conditions[-1]=conditions[-1].replace('\n','')

        
        for line in f:
            vector=line.split()

            geneName=vector[0].replace('"','')
            ratios=[float(element) for element in vector[1:]]            
            expression.append(ratios)
            
    E=numpy.array(expression)

    return E,conditions

def dataAnalyser(E,conditions):

    '''
    analysedData[condition][Q1/median/Q3]=[]
    '''

    analysedData={}

    # defining Q1s, medians and Q3s per condition
    for i in range(E.shape[1]):
        localExpression=E[:,i]
        localCondition=sampleMetadata[conditions[i]]

        if localCondition not in analysedData.keys():
            analysedData[localCondition]={}
            analysedData[localCondition]['Q1']=[]
            analysedData[localCondition]['median']=[]
            analysedData[localCondition]['Q3']=[]

        analysedData[localCondition]['Q1'].append(numpy.percentile(localExpression,25))
        analysedData[localCondition]['median'].append(numpy.median(localExpression))
        analysedData[localCondition]['Q3'].append(numpy.percentile(localExpression,75))
        
    # defining averages per condition
    broadData={}
    for condition in analysedData.keys():
        broadData[condition]=numpy.mean(analysedData[condition]['median'])
    sortedConditions=sorted(broadData,key=broadData.__getitem__)

    return analysedData,sortedConditions

def metadataReader():

    sampleMetadata={}
    annotationValues=[]
    
    with open(metadataFile,'r') as f:
        next(f)
        next(f)
        for line in f:
            vector=line.split('\t')
            sampleName=vector[0]
            annotation=vector[1]
            sampleMetadata[sampleName]=annotation

            if annotation not in annotationValues:
                annotationValues.append(annotation)

    annotationValues.sort()

    return sampleMetadata,annotationValues

def plotter(label,analysedData,sortedConditions):

    figureFile='figures/{}.pdf'.format(label)
    customMap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0,vmax=len(sortedConditions)),cmap='tab20')
    shift=0
    ceil=0; floor=0
    legendHandles=[]
    
    # selecting the sorted conditions
    for broadCondition in sortedConditions:
        
        # sorting the medians
        medians=analysedData[broadCondition]['median']
        sortedMedians=copy.deepcopy(medians)
        sortedMedians.sort()

        # obtaining the sorted quantiles
        sortedIndex=numpy.argsort(medians)
        bottom=[analysedData[broadCondition]['Q1'][index] for index in sortedIndex]
        top=[analysedData[broadCondition]['Q3'][index] for index in sortedIndex]

        # obtaining vertical ranges
        if max(top) > ceil:
            ceil=max(top)
        if min(bottom) < floor:
            floor=min(bottom)
    
        # obtaining the color and legend handles
        theColor=customMap.to_rgba(sortedConditions.index(broadCondition))
        legendHandles.append(matplotlib.patches.Patch(color=theColor,label=broadCondition))

        # plotting
        x=[i+shift for i in range(len(medians))]
        shift=x[-1]+1
        y=sortedMedians

        for i in range(len(x)):
            matplotlib.pyplot.plot([x[i],x[i]],[bottom[i],top[i]],ls='-',lw=0.355,color=theColor)
        for i in range(len(x)):
            matplotlib.pyplot.plot([x[i]],[y[i]],marker='o',mew=0,mfc='black',ms=0.4)

    # building a legend with phantom data
    matplotlib.pyplot.legend(handles=legendHandles,ncol=5,fontsize=8,columnspacing=0.2,handletextpad=0.1,loc=8,bbox_to_anchor=(0.5, 1))

    # aesthetics and closing figure
    matplotlib.pyplot.plot([0,x[-1]],[0,0],ls=':',lw=1,color='black',zorder=0)
    matplotlib.pyplot.grid(alpha=0.3)
    
    rank=x[-1]; left=-0.0*rank; right=rank+0.0*rank
    matplotlib.pyplot.xlim([left,right])

    theFloor=numpy.floor(floor)
    theCeiling=numpy.ceil(ceil)
    matplotlib.pyplot.ylim([theFloor,theCeiling])
    
    matplotlib.pyplot.xlabel('condition index')
    matplotlib.pyplot.ylabel('standardized log$_2$ relative expresion')
    matplotlib.pyplot.tight_layout(rect=(0,0,1,0.92))
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()


    sys.exit()

    return None

# 0. user define variables
dataDir='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/expressionSelectedCorems/'
metadataFile='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/metadata/array_annot.txt'

# 1. setting expression files
elements=os.listdir(dataDir)
coremPaths=[element for element in elements if '.txt' in element]
coremPaths.sort()

# 2. reading sample metadata
sampleMetadata,annotationValues=metadataReader()

# 2. iterating over the corems

for case in coremPaths:

    print(case)
    
    inputFileName=dataDir+case

    # 2.1. reading data
    E,conditions=coremReader(inputFileName)

    # 2.2. analysing data
    analysedData,sortedConditions=dataAnalyser(E,conditions)
    
    # 2.3. plotting data 
    label=case.split('.txt')[0]
    plotter(label,analysedData,sortedConditions)
