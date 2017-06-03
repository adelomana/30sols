###
### this script tests the hypothesis that protein expression can be explained by trends in RPF, i.e., TLR
###

import os,sys,numpy,random
import matplotlib,matplotlib.pyplot
import pyqt_fit,pyqt_fit.nonparam_regression,pyqt_fit.npr_methods,pyqt_fit.bootstrap
import scipy,scipy.stats

#from lmfit import Model, CompositeModel

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def backgroundCalculator(size,layer,time,backgroundOutliers):

    dist=[]
    iterations=100
    if layer == 'pt':
        source=log2proteome
    elif layer == 'mRNA':
        source=log2transcriptome
    else:
        print('error, exit')
        sys.exit()

    for i in range(iterations):
        selected=random.sample(backgroundOutliers,size)
        for element in selected:
            x=[]
            for replicate in sortedReplicateLabels:
                value=None
                try:
                    if layer == 'pt':
                        value=log2proteome['lysate'][replicate][time][element]
                    elif layer == 'mRNA':
                        value=log2transcriptome['trna'][replicate][time][element]
                    else:
                        print('error, exit')
                        sys.exit()
                except:
                    pass
                if value != None:
                    x.append(value)
            if x != []:
                average=numpy.median(x)
                dist.append(numpy.median(x))

    print(size,layer,time,len(dist))
            
    return dist

def estimation(x,y):

    est=pyqt_fit.nonparam_regression.NonParamRegression(x,y,method=pyqt_fit.npr_methods.LocalPolynomialKernel(q=0))
    est.fit()

    return est

def individualPlotter(TP,PO,TO,RO,deltaRA,RHO,TLR):

    markers=['o','s','^']
    theAlpha=0.25

    theColor='blue'
    matplotlib.pyplot.plot(TP,numpy.mean(PO,1),'-',color=theColor,label='pt',lw=0.5)
    for j in range(PO.shape[1]):
        matplotlib.pyplot.plot(TP,PO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='red'
    matplotlib.pyplot.plot(TP,numpy.mean(TO,1),'-',color=theColor,label='mRNA',lw=0.5)
    for j in range(TO.shape[1]):
        matplotlib.pyplot.plot(TP,TO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='orange'
    matplotlib.pyplot.plot(TP,numpy.mean(RO,1),'-',color=theColor,label='RBF',lw=0.5)
    for j in range(RO.shape[1]):
        matplotlib.pyplot.plot(TP,RO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

    theColor='green'
    matplotlib.pyplot.plot(TP,deltaRA,'-',color=theColor,label='RO',lw=1)

    theColor='purple'
    matplotlib.pyplot.plot(TP,RHO,'-',color=theColor,label='$\\beta$',lw=1)

    theColor='black'
    matplotlib.pyplot.plot(TP,TLR,'-',color=theColor,label='TLR',lw=2)

    if name in ribosomalProteinGenes:
        figureName='figures/ribosomalProteins/{}.pdf'.format(name)
    else:
        figureName='figures/regularProteins/{}.pdf'.format(name)

    matplotlib.pyplot.xlabel('time (h)')
    matplotlib.pyplot.ylabel('log$_2$ FC')
    matplotlib.pyplot.title(name)
    matplotlib.pyplot.legend(loc=2)
    matplotlib.pyplot.xlim([-10,45])
    
    allValues=[PO,TO,RO,deltaRA,RHO,TLR]
    topValues=[numpy.amax(element) for element in allValues]
    bottomValues=[numpy.amin(element) for element in allValues]
    if max(topValues) > 1.9:
        a=max(topValues)+0.1*max(topValues)
    else:
        a=2
    if min(bottomValues) < -1.9:
        b=min(bottomValues)+0.1*min(bottomValues)
    else:
        b=-2    
    matplotlib.pyplot.ylim([b,a])
    
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

def outliersChecker(topOutliers,bottomOutliers,backgroundOutliers):

    '''
    this function checks outliers expression distribution in the next time points
    '''

    print('top',len(topOutliers))
    print('bottom',len(bottomOutliers))

    # defining expression (mRNA and pt) of top and bottom outliers in t2 and t4
    boxData={}
    boxData['mRNA.t2']=[]
    boxData['mRNA.t4']=[]
    boxData['pt.t2']=[]
    boxData['pt.t4']=[]
    for outlier in topOutliers:
        a=[]; b=[]; c=[]; d=[]
        for replicate in sortedReplicateLabels:
            a.append(log2transcriptome['trna'][replicate]['tp2vs1'][outlier])
            b.append(log2transcriptome['trna'][replicate]['tp4vs1'][outlier])

            try:
                c.append(log2proteome['lysate'][replicate]['tp2vs1'][outlier])
                d.append(log2proteome['lysate'][replicate]['tp4vs1'][outlier])
            except:
                pass

        boxData['mRNA.t2'].append(numpy.median(a))
        boxData['mRNA.t4'].append(numpy.median(b))
        if len(c) != 0:
            boxData['pt.t2'].append(numpy.median(c))
        if len(d) != 0:
            boxData['pt.t4'].append(numpy.median(d))

    o1=boxData['mRNA.t2']
    o2=boxData['pt.t2']
    o3=boxData['mRNA.t4']
    o4=boxData['pt.t4']

    b1=backgroundCalculator(len(o1),'mRNA','tp2vs1',backgroundOutliers)
    b2=backgroundCalculator(len(o2),'pt','tp2vs1',backgroundOutliers)
    b3=backgroundCalculator(len(o3),'mRNA','tp4vs1',backgroundOutliers)
    b4=backgroundCalculator(len(o4),'pt','tp4vs1',backgroundOutliers)

    saved=b1

    results1=scipy.stats.mannwhitneyu(o1,b1)
    results2=scipy.stats.mannwhitneyu(o2,b2)
    results3=scipy.stats.mannwhitneyu(o3,b3)
    results4=scipy.stats.mannwhitneyu(o4,b4)
    print(results1,numpy.median(o1),numpy.median(b1))
    print(results2)
    print(results3)
    print(results4)
    
    matplotlib.pyplot.boxplot([o1,b1,o2,b2,o3,b3,o4,b4],labels=['mRNA.t2','B','pt.t2','B','mRNA.t4','B','pt.t4','B'],positions=[1,1.75,3,3.75,5,5.75,7,7.75])
    matplotlib.pyplot.ylabel('log$_2$ FC')
    matplotlib.pyplot.title('top outliers')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('figures/topOutliers.pdf')
    matplotlib.pyplot.clf()

    ###
    ###
    ###
    
    boxData={}
    boxData['mRNA.t2']=[]
    boxData['mRNA.t4']=[]
    boxData['pt.t2']=[]
    boxData['pt.t4']=[]
    
    for outlier in bottomOutliers:
        a=[]; b=[]; c=[]; d=[]
        for replicate in sortedReplicateLabels:
            a.append(log2transcriptome['trna'][replicate]['tp2vs1'][outlier])
            b.append(log2transcriptome['trna'][replicate]['tp4vs1'][outlier])

            try:
                c.append(log2proteome['lysate'][replicate]['tp2vs1'][outlier])
                d.append(log2proteome['lysate'][replicate]['tp4vs1'][outlier])
            except:
                pass

        boxData['mRNA.t2'].append(numpy.median(a))
        boxData['mRNA.t4'].append(numpy.median(b))
        if len(c) != 0:
            boxData['pt.t2'].append(numpy.median(c))
        if len(d) != 0:
            boxData['pt.t4'].append(numpy.median(d))

    o1=boxData['mRNA.t2']
    o2=boxData['pt.t2']
    o3=boxData['mRNA.t4']
    o4=boxData['pt.t4']

    b1=backgroundCalculator(len(o1),'mRNA','tp2vs1',backgroundOutliers)
    b2=backgroundCalculator(len(o2),'pt','tp2vs1',backgroundOutliers)
    b3=backgroundCalculator(len(o3),'mRNA','tp4vs1',backgroundOutliers)
    b4=backgroundCalculator(len(o4),'pt','tp4vs1',backgroundOutliers)

    results1=scipy.stats.mannwhitneyu(o1,b1)
    results2=scipy.stats.mannwhitneyu(o2,b2)
    results3=scipy.stats.mannwhitneyu(o3,b3)
    results4=scipy.stats.mannwhitneyu(o4,b4)
    print(results1,numpy.median(o1),numpy.median(b1))
    print(results2)
    print(results3,numpy.median(o3),numpy.median(b3))
    print(results4)
    
    matplotlib.pyplot.boxplot([o1,b1,o2,b2,o3,b3,o4,b4],labels=['mRNA.t2','B','pt.t2','B','mRNA.t4','B','pt.t4','B'],positions=[1,1.7,3,3.7,5,5.7,7,7.7])
    matplotlib.pyplot.ylabel('log$_2$ FC')
    matplotlib.pyplot.title('bottom outliers')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('figures/bottom.pdf')
    matplotlib.pyplot.clf()

    check=scipy.stats.mannwhitneyu(b1,saved)
    print('negative check',check)
    
    return None

def outliersDefiner(x,y,z):

    '''
    this function detects transcripts with more or lees ribosome occupancy than expected
    '''

    grid=numpy.linspace(numpy.log10(2),3.9,num=20)
    means=[]
    pos=[]
    devs=[]
    redIndexes=[]
    blueIndexes=[]
    topOutliers=[]; bottomOutliers=[]; backgroundOutliers=[]

    for i in range(len(grid)-1):
        a=grid[i]; b=grid[i+1]
        pos.append(numpy.mean([a,b]))
        cluster=[]
        for j in range(len(x)):
            if a<=x[j] and x[j]<=b:
                cluster.append(y[j])
        means.append(numpy.mean(cluster))
        devs.append(numpy.std(cluster))

        canada=numpy.mean(cluster)+numpy.std(cluster)*1.96
        mexico=numpy.mean(cluster)-numpy.std(cluster)*1.96

        for j in range(len(x)):
            if a<=x[j] and x[j]<=b:
                if y[j] > canada:
                    redIndexes.append(j)
                    if z[j] not in topOutliers:
                        topOutliers.append(z[j])
                elif y[j] < mexico:
                    blueIndexes.append(j)
                    if z[j] not in bottomOutliers:
                        bottomOutliers.append(z[j])
                else:
                    if z[j] not in backgroundOutliers:
                        backgroundOutliers.append(z[j])

    top=numpy.array(means)+numpy.array(devs)*1.96
    bottom=numpy.array(means)-numpy.array(devs)*1.96

    redx=[x[i] for i in range(len(x)) if i in redIndexes]
    redy=[y[i] for i in range(len(x)) if i in redIndexes]

    bluex=[x[i] for i in range(len(x)) if i in blueIndexes]
    bluey=[y[i] for i in range(len(x)) if i in blueIndexes]

    blackx=[x[i] for i in range(len(x)) if i not in blueIndexes and i not in redIndexes]
    blacky=[y[i] for i in range(len(x)) if i not in blueIndexes and i not in redIndexes]

    matplotlib.pyplot.plot(blackx,blacky,'ok',alpha=0.1,mew=0)

    matplotlib.pyplot.plot(pos,means,'-g',lw=2)
    matplotlib.pyplot.plot(pos,top,':g')
    matplotlib.pyplot.plot(pos,bottom,':g')

    matplotlib.pyplot.plot(redx,redy,'or',mew=0)
    matplotlib.pyplot.plot(bluex,bluey,'ob',mew=0)

    matplotlib.pyplot.xlim([-0.2,6])
    matplotlib.pyplot.xlim([-0.2,6])
    matplotlib.pyplot.xlabel('mRNA log$_{10}$ (TPM+1)')
    matplotlib.pyplot.ylabel('RBF log$_{10}$ (TPM+1)')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig('figures/outliersDistribution.pdf')
    matplotlib.pyplot.clf()

    return topOutliers,bottomOutliers,backgroundOutliers

def proteomicsReader():

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC
    significance={} # same as data
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in data.keys():
            data[condition]={}; significance[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}; significance[condition][replicate]={}
            
        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}; significance[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                geneName=vector[0]

                a=float(vector[2])
                b=float(vector[6])
                c=float(vector[10])

                d=float(vector[2+2])
                e=float(vector[6+2])
                f=float(vector[10+2])

                data[condition][replicate]['tp2vs1'][geneName]=a
                data[condition][replicate]['tp3vs1'][geneName]=b
                data[condition][replicate]['tp4vs1'][geneName]=c

                significance[condition][replicate]['tp2vs1'][geneName]=d
                significance[condition][replicate]['tp3vs1'][geneName]=e
                significance[condition][replicate]['tp4vs1'][geneName]=f
            
    return data,significance

def riboListReader():

    riboProteins=[]
    with open(ribosomalProteinsList,'r') as f:
        for line in f:
            vector=line.split(';')
            name=vector[1].split('=')[1]
            formattedName=name.replace('_','')
            riboProteins.append(formattedName)
    uniqueList=list(set(riboProteins))

    return uniqueList

def riboPlotter():
    
    for name in ribosomalProteinGenes:

        print(name)

        # defining protein values
        proteinObservations=[]
        for ptReplicate in sortedReplicateLabels:
            series=[0] # manually adding the first value, log2 FC = 0 for first timepoint
            for ptTimepoint in sortedRelativeTimePointLabels:
                if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys(): 
                    value=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                    series.append(value)
                else:
                    series.append(float('nan'))
            proteinObservations.append(series)
        PO=numpy.array(proteinObservations).T
        if list(numpy.mean(PO,1))[-1] > 0:
            print(PO)

        # defining transcript values
        rnaObservations=[]
        for replicate in sortedReplicateLabels:
            series=[0]
            for tp in sortedRelativeTimePointLabels:
                value=log2transcriptome['trna'][replicate][tp][name]
                series.append(value)
            rnaObservations.append(series)
        TO=numpy.array(rnaObservations).T

        # actual plotting
        markers=['o','s','^']
        theAlpha=0.25

        theColor='blue'
        matplotlib.pyplot.plot(TP,numpy.mean(PO,1),'-',color=theColor,label='pt',lw=0.5)
        for j in range(PO.shape[1]):
            matplotlib.pyplot.plot(TP,PO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

        theColor='red'
        matplotlib.pyplot.plot(TP,numpy.mean(TO,1),'-',color=theColor,label='mRNA',lw=0.5)
        for j in range(TO.shape[1]):
            matplotlib.pyplot.plot(TP,TO[:,j],marker=markers[j],color=theColor,lw=0,mew=0,alpha=theAlpha,ms=8)

        figureName='figures/completeRibosomalProteins/{}.pdf'.format(name)
        matplotlib.pyplot.ylabel('log$_2$ FC')
        matplotlib.pyplot.xlabel('time (h)')
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

        #sys.exit()

    return None

def staticAnalysis_RNAprotein():

    print('analyzing mRNA vs protein relationship in a static manner...')

    for ptReplicate in log2proteome['lysate'].keys():
        for ptTimepoint in log2proteome['lysate'][ptReplicate].keys():

            matplotlib.pyplot.plot([-7,7],[-7,7],ls='--',lw=2,color='blue')
            matplotlib.pyplot.plot([0,0],[-7,7],ls='--',lw=2,color='blue')

            for name in consistentNames:
                if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys():

                    ptRatio=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                    ptSignificance=proteomeSignificance['lysate'][ptReplicate][ptTimepoint][name]
                    mRNAratio=log2transcriptome['trna'][ptReplicate][ptTimepoint][name]

                    if ptSignificance > 0.05:
                        theColor='black'

                    else:
                        theColor='red'
                        matplotlib.pyplot.plot(mRNAratio,ptRatio,'o',alpha=0.2,mew=0,color=theColor)

                    #matplotlib.pyplot.plot(mRNAratio,ptRatio,'o',alpha=0.2,mew=0,color=theColor)

            matplotlib.pyplot.xlim([-8,8])
            matplotlib.pyplot.ylim([-8,8])

            matplotlib.pyplot.xlabel('log$_2$ FC mRNA')
            matplotlib.pyplot.ylabel('log$_2$ FC protein')

            matplotlib.pyplot.tight_layout()

            matplotlib.pyplot.savefig('figures/relation.{}.{}.png'.format(ptReplicate,ptTimepoint))
            matplotlib.pyplot.clf()

    return None

def transcriptomeRelativeConverter():

    '''
    this function assumes consistent gene names for all conditions
    '''

    log2transcriptome={}
    
    for fraction in rnaExpression.keys():
        for replicate in rnaExpression[fraction].keys():
            for timepoint in rnaExpression[fraction][replicate].keys():
                if timepoint != 'tp.1':
                    for name in rnaExpression[fraction][replicate][timepoint].keys():
                        a=rnaExpression[fraction][replicate][timepoint][name]
                        b=rnaExpression[fraction][replicate]['tp.1'][name]

                        fc=(a+1)/(b+1)
                        log2FC=numpy.log2(fc)

                        stripped=timepoint.replace('.','')
                        relativeTimePoint=stripped+'vs1'

                        if fraction not in log2transcriptome.keys():
                            log2transcriptome[fraction]={}
                        if replicate not in log2transcriptome[fraction].keys():
                            log2transcriptome[fraction][replicate]={}
                        if relativeTimePoint not in log2transcriptome[fraction][replicate].keys():
                            log2transcriptome[fraction][replicate][relativeTimePoint]={}

                        log2transcriptome[fraction][replicate][relativeTimePoint][name]=log2FC

    return log2transcriptome

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

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
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data

# 0. user defined variables

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
ribosomalProteinsList='/Volumes/omics4tb/alomana/projects/TLR/data/genome/ribosomalProteins.txt'
#transcriptomicsDataFile='/Users/adriandelomana/tmp/data/expression/expressionMatrix.kallisto.txt'
#proteomicsDataFolder='/Users/adriandelomana/tmp/data/proteomics/all/'

# 0.2. variables
timepoints=[14.3,21.5,28.8,40.8] # needs to be double checked with arjun
sortedReplicateLabels=['br1','br2','br3']
sortedRelativeTimePointLabels=['tp2vs1','tp3vs1','tp4vs1']

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression=transcriptomicsReader()
log2transcriptome=transcriptomeRelativeConverter()

# 1.2. reading protein data
log2proteome,proteomeSignificance=proteomicsReader()

# 1.3. checking consistency of transcriptome and proteome names
transcriptomeNames=[]
for fraction in log2transcriptome.keys():
    for replicate in log2transcriptome[fraction].keys():
        for timepoint in log2transcriptome[fraction][replicate].keys():
            for name in log2transcriptome[fraction][replicate][timepoint].keys():
                if name not in transcriptomeNames:
                    transcriptomeNames.append(name)

proteomeNames=[]
for fraction in log2proteome.keys():
    for replicate in log2proteome[fraction].keys():
        for timepoint in log2proteome[fraction][replicate].keys():
            for name in log2proteome[fraction][replicate][timepoint].keys():
                if name not in proteomeNames:
                    proteomeNames.append(name)

print('\t found expression quantification for {} transcripts.'.format(len(transcriptomeNames)))
print('\t found expression quantification for {} proteins.'.format(len(proteomeNames)))

# define which proteins do not have transcript equivalent
consistentNames=[]
inconsistentNames=[]
for ptName in proteomeNames:
    if ptName in transcriptomeNames:
        consistentNames.append(ptName)
    else:
        inconsistentNames.append(ptName)
consistentNames.sort()
print('\t found transcriptome info for {} proteins.'.format(len(consistentNames)))
print('\t lost {} proteins for annotation discrepancies:'.format(len(inconsistentNames)))
print('\t\t {}'.format(','.join(inconsistentNames)))
print()

# 1.4. reading ribosomal proteins list
ribosomalProteinGenes=riboListReader()
ribosomalProteinGenes.sort()

# 2. building a figure of log2 mRNA versus log2 pt
#staticAnalysis_RNAprotein()

# check how many pt are good replicates, same for transcripts. Plot it as a % of detected. use scipy.stats.norm.interval(2/3) as a rule for good replicates.

# 3. computing TLR
print('computing variables for TLR...')
TP=numpy.array(timepoints)

# 3.1. computing pt vs mRNa ratio (PMR)
rankFullProteinTrajectories=0

for name in consistentNames:

    # defining protein values
    proteinObservations=[]
    for ptReplicate in sortedReplicateLabels:
        series=[0] # manually adding the first value, log2 FC = 0 for first timepoint
        for ptTimepoint in sortedRelativeTimePointLabels:
            if name in log2proteome['lysate'][ptReplicate][ptTimepoint].keys(): 
                value=log2proteome['lysate'][ptReplicate][ptTimepoint][name]
                series.append(value)
        proteinObservations.append(series)
    PO=numpy.array(proteinObservations).T

    if PO.shape == (4,3): # 4 timepoints, 3 trajectories, full data set
        rankFullProteinTrajectories=rankFullProteinTrajectories+1

        # defining transcript values
        rnaObservations=[]
        for replicate in sortedReplicateLabels:
            series=[0]
            for tp in sortedRelativeTimePointLabels:
                value=log2transcriptome['trna'][replicate][tp][name]
                series.append(value)
            rnaObservations.append(series)
        TO=numpy.array(rnaObservations).T

        # defining RBF values
        rbfObservations=[]
        for replicate in sortedReplicateLabels:
            series=[0]
            for tp in sortedRelativeTimePointLabels:
                value=log2transcriptome['rbf'][replicate][tp][name]
                series.append(value)
            rbfObservations.append(series)
        RO=numpy.array(rbfObservations).T

        # fitting data to a model
        # https://lmfit.github.io/lmfit-py/model.html
        #mod  = CompositeModel(Model(jump), Model(gaussian), convolve)

        # computing TLR
        deltaRA=numpy.mean(RO,1)-numpy.mean(TO,1)
        RHO=numpy.mean(PO,1)-numpy.mean(TO,1)
        TLR=RHO-deltaRA
        
        # plotting
        #individualPlotter(TP,PO,TO,RO,deltaRA,RHO,TLR)

print(rankFullProteinTrajectories)

# 3.1.1. plotting all information available for ribosomal proteins
#riboPlotter()
#sys.exit()

##### scatter plot of log2 FC and tpm
allX=[]; allY=[]
allZ=[]

for name in ribosomalProteinGenes:

    print(name)

    # defining transcript values
    for tp in sortedRelativeTimePointLabels:
        fc=[]; tpm=[]
        for replicate in sortedReplicateLabels:
            a=log2transcriptome['trna'][replicate][tp][name]
            newtp='tp.{}'.format(tp[2])
            b=rnaExpression['trna'][replicate][newtp][name]
            fc.append(a)
            tpm.append(numpy.log10(b+1))
        x=numpy.mean(tpm); y=numpy.mean(fc)
    allX.append(x)
    allY.append(y)

matplotlib.pyplot.plot(allX,allY,'ok',alpha=0.5,mew=0)

matplotlib.pyplot.xlabel('mRNA log$_{10}$ (TPM+1)')
matplotlib.pyplot.ylabel('mRNA log$_{2}$ FC')

matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/scatterTPM_FC.pdf')
matplotlib.pyplot.clf()

              

#####

# 3.2. testing for distribution of affinity in first points
x=[]
y=[]
z=[]

for name in transcriptomeNames:
    transcriptValues=[]
    footprintValues=[]
    for replicate in sortedReplicateLabels:
        valueT=numpy.log10(rnaExpression['trna'][replicate]['tp.1'][name]+1)
        valueF=numpy.log10(rnaExpression['rbf'][replicate]['tp.1'][name]+1)

        transcriptValues.append(valueT)
        footprintValues.append(valueF)

    mT=numpy.median(transcriptValues)
    mF=numpy.median(footprintValues)

    if mT > 0 and mF > 0:
        x.append(mT)
        y.append(mF)
        z.append(name)

# defining outliers
topOutliers,bottomOutliers,backgroundOutliers=outliersDefiner(x,y,z)

# are outliers ribosomal genes?
itop=list(set(topOutliers) & set(ribosomalProteinGenes))
ibottom=list(set(bottomOutliers) & set(ribosomalProteinGenes))
print('intersect ribosomal genes and top outliers',itop)
print('intersect ribosomal genes and bottom outliers',ibottom)


# checking that outliers have different expression behavior in the future
outliersChecker(topOutliers,bottomOutliers,backgroundOutliers)


# 
