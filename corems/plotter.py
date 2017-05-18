import os,sys,numpy,copy
import matplotlib,matplotlib.pyplot,matplotlib.cm,matplotlib.patches

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def backgroundComputer():

    """
    this function computes the average Q1/median/Q3 of n random set of genes (not in corem), obtaining the same structure as analysedData
    """

    background={}
    
    # define the genes in the corem, complementary genes
    numberOfGenes=len(genes)
    complementaryGenes=[element for element in allGenes if element not in genes]
    
    # for each block of conditions
    for block in sortedConditions:

        background[block]={}
        background[block]['Q1']=None
        background[block]['median']=None
        background[block]['Q3']=None

        ensemble={}
        ensemble['Q1']=[]
        ensemble['median']=[]
        ensemble['Q3']=[]

        # define the conditions to operate with
        selectedConditions=[condition for condition in conditions if sampleMetadata[condition] == block]

        # define placeholders for 
        for iteration in range(iterations):

            # select a random set of complementary genes
            geneSet=numpy.random.choice(complementaryGenes,numberOfGenes)

            # defining background expression
            white=[]
            for gene in geneSet:
                v=[]
                for selectedCondition in selectedConditions:
                    value=fullExpression[selectedCondition][gene]
                    v.append(value)
                white.append(v)
            W=numpy.array(white)

            # computing an instance of statistics and appending
            q1=numpy.percentile(W,25,axis=0)
            medians=numpy.median(W,axis=0)
            q3=numpy.percentile(W,75,axis=0)

            sortedMedians=copy.deepcopy(medians)
            sortedMedians.sort()
            sortedIndex=numpy.argsort(medians)
            bottom=[q1[index] for index in sortedIndex]
            top=[q3[index] for index in sortedIndex]
            
            ensemble['Q1'].append(bottom)
            ensemble['median'].append(sortedMedians)
            ensemble['Q3'].append(top)
                        
        # computing the averages
        background[block]['Q1']=numpy.mean(numpy.array([case for case in ensemble['Q1']]),axis=0)
        background[block]['median']=numpy.mean(numpy.array([case for case in ensemble['median']]),axis=0)
        background[block]['Q3']=numpy.mean(numpy.array([case for case in ensemble['Q3']]),axis=0)

    return background

def coremReader(inputFileName):

    expression=[]
    genes=[]

    with open(inputFileName,'r') as f:
        
        header=f.readline()
        vector=header.split('\t')
        conditions=[element.replace('"','') for element in vector]
        conditions[-1]=conditions[-1].replace('\n','')

        for line in f:
            vector=line.split()

            geneName=vector[0].replace('"','')
            genes.append(geneName)
            ratios=[float(element) for element in vector[1:]]            
            expression.append(ratios)
            
    E=numpy.array(expression)

    return E,conditions,genes

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

def expressionReader():

    fullExpression={}
    allGenes=[]
    
    with open(expressionDataFile,'r') as f:
        
        header=f.readline()
        vector=header.split('\t')
        conditions=[element.replace('"','') for element in vector]
        conditions[-1]=conditions[-1].replace('\n','')

        for condition in conditions:
            fullExpression[condition]={}

        for line in f:
            vector=line.split('\t')

            geneName=vector[0].replace('"','')
            allGenes.append(geneName)
            ratios=[float(element) for element in vector[1:]]

            for i in range(len(conditions)):
                fullExpression[conditions[i]][geneName]=ratios[i]

    return fullExpression,allGenes

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

    figureFile='figures/{}.{}.pdf'.format(label,iterations)
    customMap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0,vmax=20),cmap='tab20')
    shift=0
    ceil=0; floor=0
    legendHandles=[]
    
    # selecting the sorted conditions
    for broadCondition in sortedConditions:

        # defining background values
        backgroundBottom=background[broadCondition]['Q1']
        backgroundTop=background[broadCondition]['Q3']

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
        selectedIndex=sortedConditions.index(broadCondition)
        if selectedIndex >= 14: # avoiding grey colors, indexes 14 and 15
            selectedIndex=selectedIndex+2
        theColor=customMap.to_rgba(selectedIndex)
        legendHandles.append(matplotlib.patches.Patch(color=theColor,label=broadCondition))

        # plotting
        x=[i+shift for i in range(len(medians))]
        shift=x[-1]+1

        # plotting background
        for i in range(len(x)):
            matplotlib.pyplot.plot([x[i],x[i]],[backgroundBottom[i],backgroundTop[i]],ls='-',lw=0.355,color='black',alpha=0.1)
        
        # plotting corem distributions
        for i in range(len(x)):
            matplotlib.pyplot.plot([x[i],x[i]],[bottom[i],top[i]],ls='-',lw=0.355,color=theColor)
            
        # plotting corem medians
        y=sortedMedians
        for i in range(len(x)):
            matplotlib.pyplot.plot([x[i]],[y[i]],marker='o',mew=0,mfc='black',ms=0.4)

    # building a legend with no data
    matplotlib.pyplot.legend(handles=legendHandles,ncol=6,fontsize=8,columnspacing=0.2,handletextpad=0.1,loc=8,bbox_to_anchor=(0.5, 1))

    # aesthetics and closing figure
    matplotlib.pyplot.plot([0,x[-1]],[0,0],ls=':',lw=1,color='black',zorder=0)
    matplotlib.pyplot.grid(alpha=0.3)
    
    matplotlib.pyplot.xlim([0,x[-1]])

    theFloor=numpy.floor(floor)
    theCeiling=numpy.ceil(ceil)
    matplotlib.pyplot.ylim([theFloor,theCeiling])
    matplotlib.pyplot.ylim([-3,3])
    
    matplotlib.pyplot.xlabel('condition index')
    matplotlib.pyplot.ylabel('standardized log$_2$ relative expresion')
    matplotlib.pyplot.tight_layout(rect=(0,0,1,0.92))
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()

    return None

# 0. user define variables
dataDir='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/expressionSelectedCorems/'
metadataFile='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/metadata/array_annot.txt'
expressionDataFile='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/halo_egrin2_expression_ratios.txt'
iterations=int(1e4)

# 1. setting expression files
elements=os.listdir(dataDir)
coremPaths=[element for element in elements if '.txt' in element]
coremPaths.sort()

# 2. reading sample metadata
sampleMetadata,annotationValues=metadataReader()

# 3. reading full expression data
print('reading expression data...')
fullExpression,allGenes=expressionReader()

# 4. iterating over the corems
print('working with corems...')
for case in coremPaths:

    print('\t working with corem {}...'.format(case))
    
    inputFileName=dataDir+case

    # 4.1. reading data
    print('\t reading data...')
    E,conditions,genes=coremReader(inputFileName)
    print('\t found {} genes.'.format(len(genes)))

    # 4.2. analysing data
    print('\t analyzing data...')
    analysedData,sortedConditions=dataAnalyser(E,conditions)

    # 4.3. compute background distribution
    print('\t computing background distribution...')
    background=backgroundComputer()
    
    # 4.4. plotting data
    print('\t generating figure...')
    label=case.split('.txt')[0]
    plotter(label,analysedData,sortedConditions)

    print('')

# 5. final message
print('... done.')
