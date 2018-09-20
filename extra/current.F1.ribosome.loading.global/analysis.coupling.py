###
### This script investigates the relationship between translational efficiency and other variables like expression, transcriptional regulation, half-life and transcript length.
###

### This script generates
### (Section 2) a plot of expression versus TE for each time point. Plots not used in the main manuscript nor supplementary materials. "TE.trend.tp.*.pdf"
### (Section 3) a plot of the predicted footprints learned from the model slopes. Plot used in main text. "TE.model.pdf"
### (function transcriptionalRegulationAnalyzer) a plot that assciates TE and transcriptional regulation. Plot used in main text. "TE.transcriptional.regulation.tp.*.pdf"
### (Section 4.1) a control for length of transcript and translational efficiency. "TE.control.transcript.length.pdf"
### (Section 4.2) a control for transcript half-life. "TE.control.half-life.pdf "
### (Section 5) an analysis about the relationship of expression and half-life. "expression.half-life.pdf"

import sys,math,pandas,seaborn
import matplotlib,matplotlib.pyplot
import numpy,numpy.linalg
import scipy,scipy.stats

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def DETreader():

    '''
    This function reads the DETs from sleuth and provides a dictionary as 
    DETs[timepoint][geneName]=1/-1
    DETs of abs(FC) < 2 and max expression < 10 TPMs will be excluded.
    '''

    DETs={}
    # iterate over each time point to retrieve the DETs
    for timepoint in timepoints[1:]:
        DETs[timepoint]={}            

        flag=timepoint[-1]+'1'
        fileName='sleuthResultsRNA.{}.csv'.format(flag)
        filePath=DETsDir+fileName

        with open(filePath,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                ncName=vector[1].replace('"','')
                geneName=NCsynonyms[ncName]

                # filter out abs(log2FC) < 1 or max expression < 10 TPMs
                rna0=numpy.mean([rnaExpression['trna'][replicate]['tp.1'][geneName] for replicate in replicates])
                rna1=numpy.mean([rnaExpression['trna'][replicate][timepoint][geneName] for replicate in replicates])
                log2fc=numpy.log2(rna1/rna0)

                if abs(log2fc) > 1 and numpy.max([rna0,rna1]) > 10:
                    if log2fc > 0:
                        flag=1
                    else:
                        flag=-1
                    DETs[timepoint][geneName]=flag

    return DETs

def dotColorFinder(timepoint,geneName,cloudType):

    '''
    This function defines dot colors depending on being a DET.
    '''

    if timepoint == 'tp.1':
        dotColor='noColor'
    else:
        
        if geneName not in DETs[timepoint].keys():
            dotColor='noColor'
        else:
            if DETs[timepoint][geneName] == 1:
                dotColor='red'
            elif DETs[timepoint][geneName] == -1:
                dotColor='blue'
            else:
                print('error')
                sys.exit()
    
    return dotColor

def halfLifesReader():

    '''
    This function reads the half lifes of transcripts.
    '''

    halfLifes={}
    with open(halfLifesFile,'r') as f:
        next(f)
        for line in f:
            v=line.split()
            geneName=v[0].replace('_','')
            value=float(v[1])
            halfLifes[geneName]=value

    return halfLifes

def NCsynonymsReader():

    '''
    This function builds a dictionary on the synonyms from NC to RS.
    '''

    NCsynonyms={}; transcriptLengths={}
    
    firstLine=True
    with open(transcriptomeAnnotationFile,'r') as f:
        for line in f:
            if line[0] == '>':

                # append final length
                if firstLine == False:
                    transcriptLengths[rsName]=numpy.log10(transcriptLength)

                # obtain synonyms
                vector=line.split(' ')
                ncName=vector[0].replace('>','')
                for element in vector:
                    if 'locus_tag' in element:
                        rsName=element.split('locus_tag=')[1].replace(']','').replace('_','')
                NCsynonyms[ncName]=rsName

                # reset transcript length
                transcriptLength=0
                
            else:
                firstLine=False
                v=list(line)

                transcriptLength=transcriptLength+len(v)-1
        transcriptLengths[rsName]=numpy.log10(transcriptLength)

    return NCsynonyms,transcriptLengths

def transcriptionalRegulationAnalyzer():

    '''
    This function builds a figure with the distribution of distances from expected model, across different levels of expression, for three different sets of transcripts: all, no DET, DET+ and DET-.
    '''

    print('\t\t working on transcriptional regulation association figure...')

    # f.1. build a dictionary of dictionaries with the following structure:
    # distSets[1|2|3|4|5][all|noDET|DET+|DET-]=list of y distance to regression line

    for geneName in expSet:

        # recovering information
        mRNA=expSet[geneName][0]
        r=expSet[geneName][1]

        colorLabel=expSet[geneName][2]
        if colorLabel == 'noColor':
            subset='c.noDET'
        if colorLabel == 'red':
            subset='b.DET+'
        if colorLabel == 'blue':
            subset='d.DET-'

        # deal with expression range
        ceiling=math.ceil(mRNA)
        
        # deal with distance to diagonal
        expected=m*mRNA+c
        distance=r-expected

        # append the distance into structure
        if ceiling not in distSets.keys():
            distSets[ceiling]={}

        if 'a.all' not in distSets[ceiling].keys():
            distSets[ceiling]['a.all']=[]
        if 'c.noDET' not in distSets[ceiling].keys():
            distSets[ceiling]['c.noDET']=[]
        if 'b.DET+' not in distSets[ceiling].keys():
            distSets[ceiling]['b.DET+']=[]
        if 'd.DET-' not in distSets[ceiling].keys():
            distSets[ceiling]['d.DET-']=[]

        distSets[ceiling]['a.all'].append(distance)
        distSets[ceiling][subset].append(distance)

    # f.2. build the figure
    expressionBoxes=list(distSets.keys())
    expressionBoxes.sort()
    subsets=list(distSets[expressionBoxes[0]].keys())
    subsets.sort()

    # create a dataframe for plotting with seaborn
    pos=0
    deviations=[]
    posStamps=[]
    numberEvents=[]
    for box in expressionBoxes:
        for subset in subsets:
            pos=pos+1
            # empty addition
            deviations.append(None);posStamps.append(pos)
            rank=len(distSets[box][subset])
            numberEvents.append(rank)

            # perform a test if subset equals DET-
            if subset == 'b.DET+':
                lastPlus=distSets[box][subset]
            if subset == 'c.noDET':
                lastFlat=distSets[box][subset]
            if subset == 'd.DET-':
                lastMinus=distSets[box][subset]

                if min([len(lastPlus),len(lastMinus)]) >= 10:
                    statistic,pvalue=scipy.stats.mannwhitneyu(lastPlus,lastMinus)
                    print('\t\t\t plus vs minus',len(lastPlus),len(lastMinus),statistic,pvalue)

                    statistic,pvalue=scipy.stats.mannwhitneyu(lastPlus,lastFlat)
                    print('\t\t\t plus vs noDET',len(lastPlus),len(lastFlat),statistic,pvalue)

                    statistic,pvalue=scipy.stats.mannwhitneyu(lastFlat,lastMinus)
                    print('\t\t\t noDET vs minus',len(lastFlat),len(lastMinus),statistic,pvalue)

                    print()
                    
                else:
                    print('\t\t\t bypassing test with upregulated and downregulated DET ranks',len(lastPlus),len(lastMinus))
                    print()
            
            # add values for pandas dataframe
            if rank >= 10 and box in [2,3,4]: # no boxplots for less than 10 observations
                for element in distSets[box][subset]:
                    deviations.append(element)
                    posStamps.append(pos)
        deviations.append(None)
        pos=pos+1; posStamps.append(pos)

    deviationData=list(zip(posStamps,deviations))
    df=pandas.DataFrame(data=deviationData,columns=['Subsets','Deviation'])

    # plot violin and swarm plots with seaborn
    singlePalette=['grey','red','white','blue']
    fullPalette=[];tickLabels=[]
    for i in range(expressionBoxes[-1]+1): 
        for j in range(len(singlePalette)):
            fullPalette.append(singlePalette[j])
            tickLabels.append('{}\n{}'.format(numberEvents[i*4+j],i))
        tickLabels.append('')
        fullPalette.append('green')
    
    ax=seaborn.boxplot(x='Subsets',y='Deviation',data=df,palette=fullPalette,showfliers=False)

    # aesthetics
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Subsets')
    matplotlib.pyplot.ylabel('log$_2$ Deviation')
    
    # labels
    matplotlib.pyplot.xticks([element for element in range(len(tickLabels))],tickLabels,fontsize=12)

    # text
    matplotlib.pyplot.text(9,-8.8,'$n$',ha='left',fontsize=12)
    matplotlib.pyplot.text(9,-9.4,'$\epsilon$',ha='left',fontsize=12)
    
    # this line needs to go AFTER "labels" otherwise, it's ignored... Weird.
    ax.set_xlim(9.1,23.9)
    ax.set_ylim(-8,6)

    # close figure
    figureName='figures/TE.transcriptional.regulation.{}.pdf'.format(timepoint)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    print('\t\t figure completed.')

    return None

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    geneNames=[]; timepoints=[]; replicates=[]
    
    with open(transcriptomicsDataFile,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

            if replicate not in replicates:
                replicates.append(replicate)
            if timepoint not in timepoints:
                timepoints.append(timepoint)

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
            if geneName not in geneNames:
                geneNames.append(geneName)
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data,geneNames,timepoints,replicates

###
### MAIN
###

# 0. user defined variables
theColor={}
theColor['tp.1']='red'
theColor['tp.2']='orange'
theColor['tp.3']='green'
theColor['tp.4']='blue'

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'
transcriptomeAnnotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.1.cs.NC_002608.1.fasta'
DETsDir='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/'
halfLifesFile='/Volumes/omics4tb/alomana/projects/TLR/data/halfLife/formattedHalfLifes.strains.txt'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 1.2. read DETs
NCsynonyms,transcriptLengths=NCsynonymsReader()
DETs=DETreader()

# 1.3. read half-lifes
halfLifes=halfLifesReader()

# 2. perform analysis
print('')
print('TE analysis for each time point...')

# 2.0. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. build figures on general pattern
sat=[] # needed for second part of the script, the pattern model

totalSetx=[]; totalSety=[]; totalHollowx=[]; totalHollowy=[]

# controls: transcript length and half-life
controlLength={}; controlLength['withFootprint']=[[],[]]; controlLength['withoutFootprint']=[[],[]] # need for transcript length control
controlHalfLife={}; controlHalfLife['withFootprint']=[[],[]]; controlHalfLife['withoutFootprint']=[[],[]] # need for transcript half-life

# expression vs half-life analysis
expression2halfLife={}; expression2halfLife['withFootprint']=[[],[]]; expression2halfLife['withoutFootprint']=[[],[]]

for timepoint in timepoints:

    figureName='figures/TE.trend.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    cloudColors=[]; groundColors=[]

    # necessary for transcriptional regulation analysis figure
    expSet={}
    distSets={} # to be constructed inside transcriptionalRegulationAnalyzer function
    
    for geneName in geneNames:

        transcriptLength=transcriptLengths[geneName]
        if geneName in halfLifes:
            transcriptHalfLife=halfLifes[geneName]
        else:
            transcriptHalfLife=0
        
        # check consistency of mRNA
        mRNA_TPMs=[]; footprint_TPMs=[]
        for replicate in replicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate][timepoint][geneName])
            footprint_TPMs.append(rnaExpression['rbf'][replicate][timepoint][geneName])
            
        # data transformations and quality check
        log2M=numpy.log2(numpy.array(mRNA_TPMs)+1)
        log10M=numpy.log10(numpy.array(mRNA_TPMs)+1)
        log2F=numpy.log2(numpy.array(footprint_TPMs)+1)

        # noise
        if numpy.max(log2M) > numpy.log2(10+1): # if expression is below 10 TPMs, don't consider noise
            sem=numpy.std(log2M)/numpy.sqrt(len(log2M))
            rsem_mRNA=sem/numpy.mean(log2M)
        else:
            rsem_RNA=0
            
        if numpy.max(log2F) > numpy.log2(10+1): # if expression is below 10 TPMs, don't consider noise
            sem=numpy.std(log2F)/numpy.sqrt(len(log2F))
            rsem_RF=sem/numpy.mean(log2F)
        else:
            rsem_RF=0
        
        # medians and ratio
        m=numpy.median(log10M);
        r=numpy.median(log2F)-numpy.median(log2M)

        # differenciate between trasncripts with or without footprints
        if numpy.median(footprint_TPMs) == 0:
            if rsem_mRNA < 0.3:

                hollowx.append(m); hollowy.append(r)
                totalHollowx.append(m); totalHollowy.append(r)

                # required for TE vs transcriptional regulation analysis
                dotColor=dotColorFinder(timepoint,geneName,'ground')
                groundColors.append(dotColor)
                expSet[geneName]=[m,r,dotColor]

                # for length control
                controlLength['withoutFootprint'][0].append(transcriptLength)
                controlLength['withoutFootprint'][1].append(r)

                # for half-live
                if transcriptHalfLife > 1 and transcriptHalfLife < 25:
                    controlHalfLife['withoutFootprint'][0].append(transcriptHalfLife)
                    controlHalfLife['withoutFootprint'][1].append(r)

                    # for expression vs half-life comparison
                    if m > 0:
                        expression2halfLife['withoutFootprint'][0].append(m)
                        expression2halfLife['withoutFootprint'][1].append(transcriptHalfLife)
                        
        else:
            if rsem_mRNA < 0.3 and rsem_RF < 0.3:        

                setx.append(m); sety.append(r)
                totalSetx.append(m); totalSety.append(r)

                # required for TE vs transcriptional regulation analysis
                dotColor=dotColorFinder(timepoint,geneName,'cloud')
                cloudColors.append(dotColor)
                expSet[geneName]=[m,r,dotColor]

                # for length control
                controlLength['withFootprint'][0].append(transcriptLength)
                controlLength['withFootprint'][1].append(r)

                # for half-live
                if transcriptHalfLife > 1 and transcriptHalfLife < 25:
                    controlHalfLife['withFootprint'][0].append(transcriptHalfLife)
                    controlHalfLife['withFootprint'][1].append(r)

                    # for expression vs half-life comparison
                    if m > 0:
                        expression2halfLife['withFootprint'][0].append(m)
                        expression2halfLife['withFootprint'][1].append(transcriptHalfLife)
               
    # perform regression analysis
    print('\t regression results:')
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(setx,sety)
    print('\t\t slope',slope)
    print('\t\t intercept',intercept)
    print('\t\t r_value',r_value)
    print('\t\t pvalue',p_value)
    print('\t\t std_err',std_err) 

    # compute for the model
    m=slope
    c=intercept
    expected=list(m*numpy.array(setx)+c)

    # predicted value
    predictedRatio=m*numpy.array(4)+c
    predictedValue=(2**predictedRatio)*10000
    print('\t\t predicted value:',predictedValue)

    # computed model
    satx=numpy.arange(0,10000,1)
        
    factor=m*numpy.log10(satx+1) + c + numpy.log2(satx+1)
    saty=(2**(factor))-1

    sat.append([list(satx),list(saty)])
        
    # plot figure
    for i in range(len(setx)):
        matplotlib.pyplot.plot(setx[i],sety[i],'o',alpha=0.5,mew=0,color='black')
    for i in range(len(hollowx)):
        matplotlib.pyplot.plot(hollowx[i],hollowy[i],'o',alpha=0.5,mew=0,color='tan')
    
    matplotlib.pyplot.plot(setx,expected,'-',lw=1,color=theColor[timepoint])
    
    matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
    matplotlib.pyplot.ylabel('footprint/mRNA [log$_{2}$ ratio]')

    matplotlib.pyplot.xlim([-0.1,5.3])
    matplotlib.pyplot.ylim([-15.2,8.4])

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    print('\t completed analysis per time point.\n')

    # call function to create PLOT3, i.e., TE vs transcriptional regulation
    if timepoint != 'tp.1':
        transcriptionalRegulationAnalyzer()

print('building figure for all time points...')
matplotlib.pyplot.plot(totalSetx,totalSety,'o',alpha=0.0333,mew=0,color='black')
matplotlib.pyplot.plot(totalHollowx,totalHollowy,'o',alpha=0.0333,mew=0,color='tan')

print('\t regression results:')
slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(totalSetx,totalSety)
print('\t\t slope',slope)
print('\t\t intercept',intercept)
print('\t\t r_value',r_value)
print('\t\t pvalue',p_value)
print('\t\t std_err',std_err)

# compute for the model
m=slope
c=intercept
expected=list(m*numpy.array(totalSetx)+c)

matplotlib.pyplot.plot(totalSetx,expected,'-',lw=2,color='black')

matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.trend.all.pdf')
matplotlib.pyplot.clf()
print('general model plotted.')
print('')

# 3. plot pattern from model, or panel b
print('working on infered model predictions...')

figureName='figures/TE.model.pdf'
theColors=['red','orange','green','blue']

for i in range(len(sat)):
    px=sat[i][0]
    py=sat[i][1]
    matplotlib.pyplot.plot(px,py,'-',lw=2,color=theColors[i])
    print('\t last value',py[-1])

matplotlib.pyplot.xlabel('mRNA [TPM]')
matplotlib.pyplot.ylabel('predicted footprint [TPM]')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()
print('model predictions plotted.')

# 4. controls: TE vs transcript length and half-life
print()
print('working with TE controls...')

# 4.1. TE and transcript length control
print('\t working on TE and transcript length...')

x=controlLength['withoutFootprint'][0]
y=controlLength['withoutFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='tan')

x=controlLength['withFootprint'][0]
y=controlLength['withFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='black')

slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
print('\t\t slope',slope)
print('\t\t intercept',intercept)
print('\t\t r_value',r_value)
print('\t\t pvalue',p_value)
print('\t\t std_err',std_err) 
m=slope
c=intercept
sequence=numpy.arange(min(x),max(x)+0.1,0.1)
expected=list(m*numpy.array(sequence)+c)
matplotlib.pyplot.plot(sequence,expected,'-',lw=2,color='black')

matplotlib.pyplot.xlabel('transcript length (log$_{10}$ bp)')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.control.transcript.length.pdf')
matplotlib.pyplot.clf()

# 4.2. TE and half-life control
print('\t working on TE and half-life...')
x=controlHalfLife['withoutFootprint'][0]
y=controlHalfLife['withoutFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='tan')

x=controlHalfLife['withFootprint'][0]
y=controlHalfLife['withFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='black')

slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
print('\t\t slope',slope)
print('\t\t intercept',intercept)
print('\t\t r_value',r_value)
print('\t\t pvalue',p_value)
print('\t\t std_err',std_err) 
m=slope
c=intercept
sequence=numpy.arange(min(x),max(x)+0.1,0.1)
expected=list(m*numpy.array(sequence)+c)
matplotlib.pyplot.plot(sequence,expected,'-',lw=2,color='black')

matplotlib.pyplot.xlabel('half life (min)')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.control.half-life.pdf')
matplotlib.pyplot.clf()

# 5. relationship between expression and half-life
print()
print('working on expression versus half-life relationship...')

print('\t working with transcripts without footprints...')
x=expression2halfLife['withoutFootprint'][0]
y=expression2halfLife['withoutFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.05,mew=0,color='tan')

slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
print('\t\t slope',slope)
print('\t\t intercept',intercept)
print('\t\t r_value',r_value)
print('\t\t pvalue',p_value)
print('\t\t std_err',std_err) 
m=slope
c=intercept
sequence=numpy.arange(min(x),max(x)+0.1,0.1)
expected=list(m*numpy.array(sequence)+c)
matplotlib.pyplot.plot(sequence,expected,'-',lw=2,color='tan')

print('\t working with transcripts with footprints...')
x=expression2halfLife['withFootprint'][0]
y=expression2halfLife['withFootprint'][1]
matplotlib.pyplot.plot(x,y,'o',alpha=0.05,mew=0,color='black')

slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
print('\t\t slope',slope)
print('\t\t intercept',intercept)
print('\t\t r_value',r_value)
print('\t\t pvalue',p_value)
print('\t\t std_err',std_err) 
m=slope
c=intercept
sequence=numpy.arange(min(x),max(x)+0.1,0.1)
expected=list(m*numpy.array(sequence)+c)
matplotlib.pyplot.plot(sequence,expected,'-',lw=2,color='black')

matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
matplotlib.pyplot.ylabel('half life (min)')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
matplotlib.pyplot.xlim([-0.1,4.9])
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/expression.half-life.pdf')
matplotlib.pyplot.clf()
