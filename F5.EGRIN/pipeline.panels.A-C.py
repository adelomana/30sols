import os,sys,numpy,copy,pickle
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
    for block in sortedBlockConditions:

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

def differencesAssessment(block,theColors):

    '''
    this function checks that the size of a block is appropriate
    '''

    print('\t{}...'.format(block))

    # define working conditions
    workingConditions=[]
    for specific in sampleMetadata:
        if sampleMetadata[specific] == block:
            workingConditions.append(specific)

    # compute F and G
    F=numpy.ones([len(theColors),len(theColors)]) # KS-based differences
    G=numpy.ones([len(theColors),len(theColors)]) # SC-based differences
    H=numpy.ones([len(theColors),len(theColors)]) # both KS and SC
    
    if len(workingConditions) > 50:
        F,G,H=similarityMatrixComputer(block,theColors,workingConditions,F,G,H)
                    
    return F,G,H,len(workingConditions)

def dimensionalityReductionAnalyses():

    '''
    this function computes PCA and tSNE for specified corems a PCA plot based on corem expression
    '''

    # 1. working with all corems
    theColors=[]; theAlphas=[]
    allFiles=os.listdir(allCoremsExpressionDir)
    coremLabels=[int(element.split('.txt')[0]) for element in allFiles if '.txt' in element]
    coremLabels.sort()

    for i in range(len(coremLabels)):
    
        if coremLabels[i] in clusteredCorems['green']:
            theColors.append('green'); theAlphas.append(.8)
        elif coremLabels[i] in clusteredCorems['red']:
            theColors.append('red'); theAlphas.append(.8)
        elif coremLabels[i] in clusteredCorems['magenta']:
            theColors.append('magenta'); theAlphas.append(.8)
        elif coremLabels[i] in clusteredCorems['blue']:
            theColors.append('blue'); theAlphas.append(.8)
        else:
            theColors.append('black'); theAlphas.append(0.1)
                
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
    pcaCaller(N,theColors,theAlphas,figureFile)
    
    # 1.4. t-SNE 
    figureFile='figures/figure.tSNE.pdf'
    perplexityValue=50
    print(N.shape,len(coremLabels))
    tSNEcaller(N,theColors,theAlphas,figureFile,perplexityValue)

    """
    
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

    """
    
    return None

def exhaustiveDifferencesFinder():

    '''
    this function computes the frequency of differences between corems of different clusters (colors)
    '''

    theColors=['red','green','blue','magenta']
    similarityJar=jarDir+'similarity.pickle'

    # compute and store similarity based on KS
    similarityKS={}
    similaritySC={}
    similarityBoth={}
    conditionsSizes={}
    for block in sortedBlockConditions:
        # compute the matrix of differences for that particular condition
        F,G,H,conditionSize=differencesAssessment(block,theColors)
        if numpy.mean(F) != 1:
            similarityKS[block]=F
            similaritySC[block]=G
            similarityBoth[block]=H
            conditionsSizes[block]=conditionSize
    similarities=[similarityKS,similaritySC,similarityBoth,conditionsSizes]
    f=open(similarityJar,'wb')
    pickle.dump(similarities,f)
    f.close()
    
    # recover similarity
    f=open(similarityJar,'rb')
    similarities=pickle.load(f)
    f.close()

    # plot a heat map with similarities
    heatmapSimilaritiesPlotter(similarities)
    
    return None

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

def geneMembershipReader():

    '''
    this function reads the gene members of selected corems
    '''

    coremGeneMemberships={}
    allFiles=os.listdir(allCoremsExpressionDir)
    coremLabels=[int(element.split('.txt')[0]) for element in allFiles if '.txt' in element]
    coremLabels.sort()
        
    for i in range(len(coremLabels)):
        fileName='{}{}.txt'.format(allCoremsExpressionDir,coremLabels[i])
        genes=[]
        with open(fileName,'r') as f:
            next(f)
            for line in f:
                vector=line.split('\t')
                geneName=vector[0].replace('"','')
                genes.append(geneName)
        coremGeneMemberships[coremLabels[i]]=genes
            
    return coremGeneMemberships

def heatmapSimilaritiesPlotter(similarities):

    '''
    this function plots a heatmap of similarity indexes, both for Kolmogorov-Smirnov-based and Spearman Coefficient-based similarities
    '''

    methodLabels=['ks','sc','both']
    computedBlocks=list(similarities[0].keys())
    conditionsSizes=similarities[-1]

    for i in range(len(methodLabels)):

        # sorting conditions from high similarity to low
        similarityOrder={}
        
        for block in computedBlocks:
            B=similarities[i][block]
            similarityOrder[block]=numpy.mean(B)

        sortedComputedBlockConditions=sorted(similarityOrder,key=similarityOrder.__getitem__,reverse=True)

        for j in range(len(sortedComputedBlockConditions)):
            if j == 0:
                V=similarities[i][sortedComputedBlockConditions[j]]
            else:
                V=numpy.hstack((V,similarities[i][sortedComputedBlockConditions[j]]))

        # plotting figure
        figureName='figures/similarity.{}.png'.format(methodLabels[i])
    
        matplotlib.pyplot.imshow(V,interpolation='none',cmap='viridis',vmin=0.,vmax=1.)

        cb=matplotlib.pyplot.colorbar(label='similarity',orientation='horizontal',fraction=0.025) 
        cb.ax.tick_params(labelsize=10)

        positions=[(i*4)-0.5 for i in range(len(sortedComputedBlockConditions))]
        values=[str(int(element+0.5)) for element in positions]
        matplotlib.pyplot.xticks(positions,values,size=10)
        matplotlib.pyplot.yticks([],[])

        for j in range(len(sortedComputedBlockConditions)):
            xpos=2+j*4
            ypos=-1
            
            label='{} n={}'.format(sortedComputedBlockConditions[j],conditionsSizes[sortedComputedBlockConditions[j]])
            matplotlib.pyplot.text(xpos,ypos,label,rotation=90,horizontalalignment='center',verticalalignment='bottom')
    
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

    return None

def metadataReader():

    '''
    this function reads the metadata
    '''

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

def pcaCaller(N,theColors,theAlphas,figureFile):

    '''
    this function makes a PCA and its figure
    '''
    
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

def plotter(label,analysedData,sortedBlockConditions):

    figureFile='figures/{}.iter{}.pdf'.format(label,iterations)
    customMap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0,vmax=20),cmap='tab20')
    shift=0
    ceil=0; floor=0
    legendHandles=[]
    
    # selecting the sorted conditions
    for block in sortedBlockConditions:

        # defining background values
        backgroundBottom=background[block]['Q1']
        backgroundTop=background[block]['Q3']

        # sorting the medians
        medians=analysedData[block]['median']
        sortedMedians=copy.deepcopy(medians)
        sortedMedians.sort()

        # obtaining the sorted quantiles
        sortedIndex=numpy.argsort(medians)
        bottom=[analysedData[block]['Q1'][index] for index in sortedIndex]
        top=[analysedData[block]['Q3'][index] for index in sortedIndex]

        # defining the same order of conditions for background
        sortedBackgroundBottom=[backgroundBottom[index] for index in sortedIndex]
        sortedBackgroundTop=[backgroundTop[index] for index in sortedIndex]

        # obtaining vertical ranges
        if max(top) > ceil:
            ceil=max(top)
        if min(bottom) < floor:
            floor=min(bottom)
    
        # obtaining the color and legend handles
        selectedIndex=sortedBlockConditions.index(block)
        if selectedIndex >= 14: # avoiding grey colors, indexes 14 and 15
            selectedIndex=selectedIndex+2
        theColor=customMap.to_rgba(selectedIndex)
        #legendHandles.append(matplotlib.patches.Patch(color=theColor,label=broadCondition))

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

def similarityMatrixComputer(block,theColors,workingConditions,F,G,H):

    '''
    this function computes a matrix of differences for clusters of corems for a particular condition
    '''

    for i in range(len(theColors)):
        for j in range(len(theColors)):
            if i<=j:
                color1=theColors[i]
                color2=theColors[j]
                
                # define corems
                coremsA=clusteredCorems[color1]
                coremsB=clusteredCorems[color2]

                # define frequency of differences
                maxRankDifferences=0
                KSdifferences=0
                SCdifferences=0
                anyDifferences=0
                
                for coremA in coremsA:
                    for coremB in coremsB:

                        # obtain median expression for coremA
                        X=[]
                        genesA=coremGeneMemberships[coremA]
                        for gene in genesA:
                            Y=[]
                            for condition in workingConditions:
                                Y.append(fullExpression[condition][gene])
                            X.append(Y)
                        Z=numpy.array(X)
                        medianA=numpy.median(Z,axis=0)

                        # obtain median expression for coremB
                        X=[]
                        genesB=coremGeneMemberships[coremB]
                        for gene in genesB:
                            Y=[]
                            for condition in workingConditions:
                                Y.append(fullExpression[condition][gene])
                            X.append(Y)
                        Z=numpy.array(X)
                        medianB=numpy.median(Z,axis=0)

                        # test for the difference between them: KS and Spearman rank
                        maxRankDifferences=maxRankDifferences+1
                        diffKS=0
                        diffSC=0

                        D,pvalueKS=scipy.stats.ks_2samp(medianA,medianB)
                        if pvalueKS <= 0.05:
                            KSdifferences=KSdifferences+1
                            diffKS=1

                        correlation,pvalueSC=scipy.stats.spearmanr(medianA,medianB)
                        
                        if correlation >= 0.4 and pvalueSC <= 0.05: 
                            pass
                        else:
                            SCdifferences=SCdifferences+1
                            diffSC=1

                        if diffKS == 1 or diffSC == 1:
                            anyDifferences=anyDifferences+1


                        '''
                        
                        # plotting rank figures
                        if diffKS == 0 and diffSC == 1:
                            figureName='figures/similarityFigures/{}.{}.{}.{}.{}.pdf'.format(block,color1,color2,coremA,coremB)

                            A=numpy.vstack([medianA,numpy.ones(len(medianA))]).T
                            m,c=numpy.linalg.lstsq(A,numpy.array(medianB))[0]

                            matplotlib.pyplot.scatter(medianA,medianB,color='black')
                            matplotlib.pyplot.plot(medianA,medianA*m+c,color='red',lw=2)

                            message='D={:.2f}, p={:.3f}, diff={:d}\nSC={:.2f}, p={:.3f}, diff={:d}'.format(D,pvalueKS,diffKS,correlation,pvalueSC,diffSC)
                            matplotlib.pyplot.text(0,1,message)
                            
                            print('***\t {} {}'.format(coremA,coremB))
                            print(message)
                            print(figureName)
                            print('')

                            matplotlib.pyplot.tight_layout()
                            matplotlib.pyplot.savefig(figureName)
                            matplotlib.pyplot.clf()
                        '''

                        
                # define similarity
                similarityKS=1-(KSdifferences/maxRankDifferences)
                F[i,j]=similarityKS
                F[j,i]=similarityKS

                similaritySC=1-(SCdifferences/maxRankDifferences)
                G[i,j]=similaritySC
                G[j,i]=similaritySC
                
                similarityBoth=1-(anyDifferences/maxRankDifferences)
                H[i,j]=similarityBoth
                H[j,i]=similarityBoth

                print('\t\t {} {} {}'.format(color1,color2,similarityBoth))

    return F,G,H

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
allCoremsExpressionDir='/Volumes/omics4tb/alomana/projects/TLR/data/HsaEGRIN/allCorems/'
dataDir='/Volumes/omics4tb/alomana/projects/TLR/data/HsaEGRIN/expressionSelectedCorems/'
metadataFile='/Volumes/omics4tb/alomana/projects/TLR/data/HsaEGRIN/metadata/array_annot.txt'
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/HsaEGRIN/halo_egrin2_expression_ratios.txt'
jarDir='/Volumes/omics4tb/alomana/projects/TLR/data/HsaEGRIN/jars/'
iterations=int(1e1)

clusteredCorems={}
clusteredCorems['green']=[8655,12578,7474,7473,20960,20806,7469]
clusteredCorems['red']=[14306,2822,9692,21843,7884,4320,29824,21846,3882,17167,7870,3570,4317,472,2383,2887,20778,3415,8640,7321,1738]
clusteredCorems['magenta']=[1842,8639,4277,4084,3418,1852,1845,1841,1837,1674,1808,2978,8630,753,7317,20780,3403,2980,3597,3428,3278,31951,2976,21855,0,1143]
clusteredCorems['blue']=[21857,3404,7316,20779,639,3280,164,31942]

# 1. reading data
# 1.1. setting expression files
elements=os.listdir(dataDir)
coremPaths=[element for element in elements if '.txt' in element]
coremPaths.sort()

# 1.2. reading sample metadata
sampleMetadata,annotationValues=metadataReader()

# 1.3. reading full expression data
print('reading expression data...')
fullExpression,allGenes,allConditions=expressionReader()

# 1.4. reading corem gene memberships for clustered ones
coremGeneMemberships=geneMembershipReader()

# 2. iterating over specific corems
print('working with corems...')
aggregateData={}
for case in coremPaths:

    print('\t working with corem {}...'.format(case))
    
    inputFileName=dataDir+case

    # 2.1. reading data
    print('\t reading data...')
    E,conditions,genes=coremReader(inputFileName)
    print('\t found {} genes.'.format(len(genes)))
    aggregateData[case]=(E,conditions)

    # 2.2. analysing data
    print('\t analyzing data...')
    analysedData,sortedBlockConditions=dataAnalyser(E,conditions)

    # 2.3. compute background distribution
    print('\t computing background distribution...')
    background=backgroundComputer()
    
    # 2.4. plotting data
    print('\t generating figure...')
    label=case.split('.txt')[0]
    plotter(label,analysedData,sortedBlockConditions)

    print('')

# 3. working with corem median profiles
coremAverageTrends(aggregateData)

# 4. dimensionality reduction analyses
dimensionalityReductionAnalyses()

# 5. condition-specific finder
print('computing differences between clusters of corems at specific conditions...')
exhaustiveDifferencesFinder()

# 6. final message
print('... done.')
