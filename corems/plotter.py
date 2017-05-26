import os,sys,numpy,copy
import matplotlib,matplotlib.pyplot,matplotlib.cm,matplotlib.patches
import scipy,scipy.stats
import sklearn,sklearn.decomposition,sklearn.manifold

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

            ensemble['Q1'].append(q1)
            ensemble['median'].append(medians)
            ensemble['Q3'].append(q3)
                        
        # computing the averages
        background[block]['Q1']=numpy.mean(numpy.array([case for case in ensemble['Q1']]),axis=0)
        background[block]['median']=numpy.mean(numpy.array([case for case in ensemble['median']]),axis=0)
        background[block]['Q3']=numpy.mean(numpy.array([case for case in ensemble['Q3']]),axis=0)

    return background

def coremAverageTrends(aggregateData):

    '''
    this function computes the median expression of common conditions for corems of interest.
    then it computes the percentile-based accumulative distribution to assess patterns of expression differences
    '''

    ### f.1. data organization
    # f.1.1. define corem names
    coremNames=list(aggregateData.keys())
    coremNames.sort()

    # f.1.2. define the common conditions in dictionary
    allConditions=[set(aggregateData[name][1]) for name in coremNames]
    commons=list(set.intersection(*allConditions))

    commonConditions={}
    for common in commons:
        broad=sampleMetadata[common]
        if broad not in commonConditions.keys():
            commonConditions[broad]=[common]
        else:
            commonConditions[broad].append(common)

    broadConditions=list(commonConditions.keys())
    sortedBroadConditions=sorted(broadConditions,key=str.lower)

    # f.1.3. obtain median data for common conditions sorting and computing cumulative distribution
    Z=[]
    for name in coremNames:
        W=[]
        for condition in commons:
            conditionIndex=aggregateData[name][1].index(condition)
            medianValue=numpy.median(aggregateData[name][0][:,conditionIndex])
            W.append(medianValue)
        Z.append(W)

    # f.1.4. sorting expression based on first corem
    sortedIndexes=numpy.argsort(Z[0])

    ### f.2. plots
    x=[i for i in range(len(sortedIndexes))]
    # ['1738.txt', '1842.txt', '639.txt', '7469.txt'] should be gene module indexes
    # 72, 10, 43, 9 should be colors
    # red, magenta, blue, green
    theColors=['red','magenta','blue','green']
    
    # f.2.1. expression plots
    figureFile='figures/corem.expression.pdf'
    for i in range(len(Z)):
        sortedElements=[Z[i][index] for index in sortedIndexes]
        matplotlib.pyplot.plot(x,sortedElements,'-',color=theColors[i],alpha=0.5,label=theColors[i],lw=1,zorder=len(Z)-i)
    matplotlib.pyplot.xlabel('condition')
    matplotlib.pyplot.ylabel('standardized log$_2$ relative expresion')
    matplotlib.pyplot.legend(loc=4,ncol=2,fontsize=12)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()
    
    # f.2.2. percentile plots
    figureFile='figures/corem.percentiles.pdf'
    for i in range(len(Z)):
        sortedElements=[Z[i][index] for index in sortedIndexes]
        percentiles=(numpy.argsort(sortedElements)/len(sortedElements))*100
        matplotlib.pyplot.plot(x,percentiles,'.',color=theColors[i],alpha=0.5,mew=0.,ms=15,label=theColors[i],zorder=len(Z)-i)
    matplotlib.pyplot.xlabel('condition')
    matplotlib.pyplot.ylabel('condition percentile')
    matplotlib.pyplot.legend(loc=4,ncol=2,fontsize=12)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()

    # f.2.3. median expression plot per condition
    figureFile='figures/coremMediaProfiles.pdf'
    shift=0
    
    for broad in sortedBroadConditions:
        print(broad)
        T=[]
        for name in coremNames:
            S=[]
            for specific in commonConditions[broad]:
                conditionIndex=aggregateData[name][1].index(specific)
                medianValue=numpy.median(aggregateData[name][0][:,conditionIndex])
                S.append(medianValue)
            #W.sort()
            T.append(S)
        
        for i in range(len(T)):
            x=[i+shift for i in range(len(T[i]))]
            sortedIndexes=numpy.argsort(T[0])
            y=[T[i][index] for index in sortedIndexes]
            matplotlib.pyplot.plot(x,y,'-',color=theColors[i],lw=0.5)
        
        shift=shift+len(T[i])

        # calling Kruskal
        (statistic,pvalue)=scipy.stats.kruskal(T[0],T[1],T[2],T[3])
        print(statistic,pvalue)
        if pvalue < 0.05:
             matplotlib.pyplot.plot([shift],[2.5],'*k')
            
    matplotlib.pyplot.ylim([-4.3,3])
    matplotlib.pyplot.plot([0,shift],[0,0],ls=':',lw=1,color='black',zorder=0)
    matplotlib.pyplot.grid(alpha=0.3)
    matplotlib.pyplot.xlabel('condition index')
    matplotlib.pyplot.ylabel('standardized log$_2$ relative expresion')
    matplotlib.pyplot.tight_layout(rect=(0,0,1,0.92))
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()



    # f.2.4. PCA

    # correlation plots. NB asks for background. compare with noise the rank correlation

    return None

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
        
    # defining conditions order alphabetically
    conditions=list(analysedData.keys())
    sortedConditions=sorted(conditions,key=str.lower)

    return analysedData,sortedConditions

def expressionReader():

    '''
    this function creates a dictionary with all expression values as expression[condition][geneName]=value
    '''

    fullExpression={}
    allGenes=[]
    
    with open(expressionDataFile,'r') as f:
        
        header=f.readline()
        vector=header.split('\t')
        allConditions=[element.replace('"','') for element in vector]
        allConditions[-1]=allConditions[-1].replace('\n','')

        for condition in allConditions:
            fullExpression[condition]={}

        for line in f:
            vector=line.split('\t')

            geneName=vector[0].replace('"','')
            allGenes.append(geneName)
            ratios=[float(element) for element in vector[1:]]

            for i in range(len(allConditions)):
                fullExpression[allConditions[i]][geneName]=ratios[i]

    return fullExpression,allGenes,allConditions

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

def dimensionalityReductionAnalyses():

    '''
    this function computes PCA and tSNE for specified corems a PCA plot based on corem expression
    '''

    # 1. working with all corems
    theColors=[]; theAlphas=[]
    # 1.1. reading gene memberships for corems
    coremGeneMemberships={}
    allFiles=os.listdir(allCoremsExpressionDir)
    coremLabels=[int(element.split('.txt')[0]) for element in allFiles if '.txt' in element]
    coremLabels.sort()
        
    for i in range(len(coremLabels)):
        
        if coremLabels[i] in greenLabels:
            theColors.append('green'); theAlphas.append(.8)
        elif coremLabels[i] in redLabels:
            theColors.append('red'); theAlphas.append(.8)
        elif coremLabels[i] in magentaLabels:
            theColors.append('magenta'); theAlphas.append(.8)
        elif coremLabels[i] in blueLabels:
            theColors.append('blue'); theAlphas.append(.8)
        else:
            theColors.append('black'); theAlphas.append(0.1)
            
        fileName='{}{}.txt'.format(allCoremsExpressionDir,coremLabels[i])
        genes=[]
        with open(fileName,'r') as f:
            next(f)
            for line in f:
                vector=line.split('\t')
                geneName=vector[0].replace('"','')
                genes.append(geneName)
        coremGeneMemberships[coremLabels[i]]=genes
                
    # 1.2. define median expression over all conditions for each corem
    M=[] # a matrix containing the medians of corems
    for i in range(len(coremLabels)):
        X=[]
        for gene in coremGeneMemberships[coremLabels[i]]:
            Y=[]
            for condition in allConditions:
                Y.append(fullExpression[condition][gene])
            X.append(Y)
        Z=numpy.array(X)
        median=numpy.median(Z,axis=0)
        M.append(median)
    N=numpy.array(M)

    # 1.3. PCA
    figureFile='figures/figure.pca.pdf'
    perplexityValue=50
    pcaCaller(N,theColors,theAlphas,figureFile)
    
    # 1.4. t-SNE 
    figureFile='figures/figure.tSNE.pdf'
    perplexityValue=50
    tSNEcaller(N,theColors,theAlphas,figureFile,perplexityValue)


    # 2. only ribosomal corems
    theColors=[]; theAlphas=[]
    for i in range(len(allRiboCorems)):
        
        if allRiboCorems[i] in greenLabels:
            theColors.append('green'); theAlphas.append(.8)
        elif allRiboCorems[i] in redLabels:
            theColors.append('red'); theAlphas.append(.8)
        elif allRiboCorems[i] in magentaLabels:
            theColors.append('magenta'); theAlphas.append(.8)
        elif allRiboCorems[i] in blueLabels:
            theColors.append('blue'); theAlphas.append(.8)
        else:
            theColors.append('black'); theAlphas.append(0.1)
            
    M=[] # a matrix containing the medians of corems
    for i in range(len(allRiboCorems)):
        X=[]
        for gene in coremGeneMemberships[allRiboCorems[i]]:
            Y=[]
            for condition in allConditions:
                Y.append(fullExpression[condition][gene])
            X.append(Y)
        Z=numpy.array(X)
        median=numpy.median(Z,axis=0)
        M.append(median)
    N=numpy.array(M)

    # 2.1. PCA
    figureFile='figures/figure.pca.ribo.pdf'
    pcaCaller(N,theColors,theAlphas,figureFile)
    
    # 2.2. t-SNE 
    figureFile='figures/figure.tSNE.ribo.pdf'
    perplexityValue=5
    tSNEcaller(N,theColors,theAlphas,figureFile,perplexityValue)
    
    return None

def pcaCaller(N,theColors,theAlphas,figureFile):
    
    print('running PCA...')
    pcaMethod=sklearn.decomposition.PCA(n_components=5)
    pcaObject=pcaMethod.fit(N)
    new=pcaObject.transform(N)
    explainedVar=pcaObject.explained_variance_ratio_
    print('cumsum explained variance...')
    print(numpy.cumsum(explainedVar))

    for i in range(len(new)):
        matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theColors[i],alpha=theAlphas[i],s=60,lw=0)

    matplotlib.pyplot.xlabel('PCA 1 ({0:.2f} var)'.format(explainedVar[0]))
    matplotlib.pyplot.ylabel('PCA 2 ({0:.2f} var)'.format(explainedVar[1]))
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()
    print()

    return None

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

        # defining the same order of conditions for background
        sortedBackgroundBottom=[backgroundBottom[index] for index in sortedIndex]
        sortedBackgroundTop=[backgroundTop[index] for index in sortedIndex]

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
            matplotlib.pyplot.plot([x[i],x[i]],[sortedBackgroundBottom[i],sortedBackgroundTop[i]],ls='-',lw=0.355,color='black',alpha=0.1)
        
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

def tSNEcaller(N,theColors,theAlphas,figureFile,perplexityValue):

    print('running t-SNE...')
    tSNE_Method=sklearn.manifold.TSNE(method='exact',verbose=1,init='pca',perplexity=perplexityValue)
    tSNE_Object=tSNE_Method.fit(N)
    new=tSNE_Object.fit_transform(N)

    for i in range(len(new)):
        matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theColors[i],alpha=theAlphas[i],s=60,lw=0)
    matplotlib.pyplot.xlabel('tSNE 1')
    matplotlib.pyplot.ylabel('tSNE 2')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()
    print()

    return None

# 0. user define variables
allCoremsExpressionDir='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/allCorems/'
dataDir='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/expressionSelectedCorems/'
metadataFile='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/metadata/array_annot.txt'
expressionDataFile='/Users/alomana/gDrive2/projects/TLR/data/HaloEGRIN/halo_egrin2_expression_ratios.txt'
iterations=int(1e1)

greenLabels=[8655,12578,7474,7473,20960,20806,7469]
redLabels=[14306,2822,9692,21843,7884,4320,29824,21846,3882,17167,7870,3570,4317,472,2383,2887,20778,3415,8640,7321,1738]
magentaLabels=[1842,8639,4277,4084,3418,1852,1845,1841,1837,1674,1808,2978,8630,753,7317,20780,3403,2980,3597,3428,3278,31951,2976,21855,0,1143]
blueLabels=[21857,3404,7316,20779,639,3280,164,31942]
allRiboCorems=greenLabels+redLabels+magentaLabels+blueLabels

# 1. setting expression files
elements=os.listdir(dataDir)
coremPaths=[element for element in elements if '.txt' in element]
coremPaths.sort()

# 2. reading sample metadata
sampleMetadata,annotationValues=metadataReader()

# 3. reading full expression data
print('reading expression data...')
fullExpression,allGenes,allConditions=expressionReader()

# 4. iterating over the corems
print('working with corems...')
aggregateData={}
for case in coremPaths:

    print('\t working with corem {}...'.format(case))
    
    inputFileName=dataDir+case

    # 4.1. reading data
    print('\t reading data...')
    E,conditions,genes=coremReader(inputFileName)
    print('\t found {} genes.'.format(len(genes)))
    aggregateData[case]=(E,conditions)

    # 4.2. analysing data
    #print('\t analyzing data...')
    #analysedData,sortedConditions=dataAnalyser(E,conditions)

    # 4.3. compute background distribution
    #print('\t computing background distribution...')
    #background=backgroundComputer()
    
    # 4.4. plotting data
    #print('\t generating figure...')
    #label=case.split('.txt')[0]
    #plotter(label,analysedData,sortedConditions)

    #print('')

# 5. working with corem median profiles
#coremAverageTrends(aggregateData)

# 6. dimensionality reduction analyses
dimensionalityReductionAnalyses()

# 7. condition-specific finder
#exhaustiveDifferencesFinder()

# 6. final message
print('... done.')
