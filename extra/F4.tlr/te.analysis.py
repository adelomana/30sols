###
### This script generates
### (PLOT1) a plot with all time point  observations plus all slopes.
### (PLOT2) a plot with the model from the slopes.
### (PLOT3) a colored plot with differential expression (red/blue).
### (PLOT4) a control for length of transcript and translational efficiency.
### (PLOT5) a control for transcript half-life
### 

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

def transcriptSetsDistributionAnalyzer():

    '''
    This function builds a figure with the distribution of distances from expected model, across different levels of expression, for three different sets of transcripts: all, no DET, DET+ and DET-.
    '''

    print('\t working on transcript distance distributions from expected...')

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
                    print('\t\t plus vs minus',len(lastPlus),len(lastMinus),statistic,pvalue)

                    statistic,pvalue=scipy.stats.mannwhitneyu(lastPlus,lastFlat)
                    print('\t\t plus vs noDET',len(lastPlus),len(lastFlat),statistic,pvalue)

                    statistic,pvalue=scipy.stats.mannwhitneyu(lastFlat,lastMinus)
                    print('\t\t noDET vs minus',len(lastFlat),len(lastMinus),statistic,pvalue)

                    print()
                    
                else:
                    print('\t\t',len(lastPlus),len(lastMinus))
            
            # add values for pandas dataframe
            if rank >= 10 and box in [2,3,4]: # no boxplots for less than 10 observations
            #if rank >= 0: # no boxplots for less than 10 observations 
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
    
    #ax=seaborn.violinplot(x='Subsets',y='Deviation',data=df,inner=None,palette=fullPalette)
    #ax=seaborn.swarmplot(x='Subsets',y='Deviation',data=df,palette=fullPalette)
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

    # close figure
    figureName='figures/TE.blocks.{}.seaborn.pdf'.format(timepoint)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

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
print('computing the analysis...')

# 2.0. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. build figures on general pattern
sat=[] # needed for second part of the script, the pattern model

totalSetx=[]; totalSety=[]; totalHollowx=[]; totalHollowy=[]; regressions=[]
lengthControl=[]; lengthControlHollow=[]; lengthControlRegressions=[]
degControl=[]; degControlHollow=[]; degControlRegressions=[]

degExpControl=[]; degExpRegressions=[]

highestExpression=0

for timepoint in timepoints:

    figureName='figures/TE.generalTrend.color.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    lengthx=[]
    cloudColors=[]; groundColors=[]
    expSet={}
    distSets={} # to be constructed in transcriptSetsDistributionAnalyzer
    degx=[]; degy=[]
    degExpx=[]; degExpy=[]
    
    for geneName in geneNames:
        
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

        # filter useful data
        if numpy.median(footprint_TPMs) == 0:
            if rsem_mRNA < 0.3:

                hollowx.append(m); hollowy.append(r)
                totalHollowx.append(m); totalHollowy.append(r)
                lengthControlHollow.append([transcriptLengths[geneName],r])
                try:
                    degControlHollow.append([halfLifes[geneName],r])
                    degExpControl.append([halfLifes[geneName],m])
                    if halfLifes[geneName] < 25:
                        degExpx.append(m); degExpy.append(halfLifes[geneName])
                except:
                    pass
                
                dotColor=dotColorFinder(timepoint,geneName,'ground')
                groundColors.append(dotColor)
                expSet[geneName]=[m,r,dotColor]
        else:
            if rsem_mRNA < 0.3 and rsem_RF < 0.3:        

                setx.append(m); sety.append(r)
                totalSetx.append(m); totalSety.append(r)                

                lengthx.append(transcriptLengths[geneName])
                lengthControl.append([transcriptLengths[geneName],r])
                try:
                    degControl.append([halfLifes[geneName],r])
                    degx.append(halfLifes[geneName]); degy.append(r)
                    degExpControl.append([halfLifes[geneName],m])
                    if halfLifes[geneName] < 25:
                        degExpx.append(m); degExpy.append(halfLifes[geneName])
                except:
                    pass
                
                dotColor=dotColorFinder(timepoint,geneName,'cloud')
                cloudColors.append(dotColor)
                expSet[geneName]=[m,r,dotColor]

    if max(setx) > highestExpression:
        highestExpression=max(setx)
    # perform regression analysis
    print('\t regression results:')
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(setx,sety)
    print('\t\t slope',slope)
    print('\t\t intercept',intercept)
    print('\t\t r_value',r_value)
    print('\t\t pvalue',p_value)
    print('\t\t std_err',std_err)

    print('\t regression for length control:')
    slopeL,interceptL,r_valueL,p_valueL,std_errL=scipy.stats.linregress(lengthx,sety)
    print('\t\t slope',slopeL)
    print('\t\t intercept',interceptL)
    print('\t\t r_value',r_valueL)
    print('\t\t pvalue',p_valueL)
    print('\t\t std_err',std_errL)

    print('\t regression for half-life control:')
    slopeHL,interceptHL,r_valueHL,p_valueHL,std_errHL=scipy.stats.linregress(degx,degy)
    print('\t\t slope',slopeHL)
    print('\t\t intercept',interceptHL)
    print('\t\t r_value',r_valueHL)
    print('\t\t pvalue',p_valueHL)
    print('\t\t std_err',std_errHL)

    print('\t regression for half-life / expression...')
    slopeHE,interceptHE,r_valueHE,p_valueHE,std_errHE=scipy.stats.linregress(degExpx,degExpy)
    print('\t\t slope',slopeHE)
    print('\t\t intercept',interceptHE)
    print('\t\t r_value',r_valueHE)
    print('\t\t pvalue',p_valueHE)
    print('\t\t std_err',std_errHE)

    # compute for the model
    m=slope
    c=intercept
    expected=list(m*numpy.array(setx)+c)
    regressions.append([slope,intercept])

    # length model
    lengthControlRegressions.append([slopeL,interceptL])

    # degradation model
    degControlRegressions.append([slopeHL,interceptHL])

    # degradation expression model
    degExpRegressions.append([slopeHE,interceptHE])

    # computed from Matt Wall on log2
    ### satx=2**(numpy.array(setx))-1
    ### saty=(2**c)*((satx+1)**(1+m))-1

    # computed ALO on log10 for x #!!! check model
    satx=(10**numpy.array(setx))-1
    factor=m*numpy.log10(satx+1) + c + numpy.log2(satx+1)
    saty=(2**(factor))-1

    sat.append([list(satx),list(saty)])
        
    # plot figure
    for i in range(len(setx)):
        if cloudColors[i] != 'noColor':
            matplotlib.pyplot.plot(setx[i],sety[i],'o',alpha=0.5,mew=0,color=cloudColors[i])
    for i in range(len(hollowx)):
        if groundColors[i] != 'noColor':
            matplotlib.pyplot.plot(hollowx[i],hollowy[i],'o',alpha=0.5,mew=0,color=groundColors[i])
    
    matplotlib.pyplot.plot(setx,expected,'-',lw=1,color=theColor[timepoint])
    
    matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
    matplotlib.pyplot.ylabel('footprint/mRNA [log$_{2}$ ratio]')

    matplotlib.pyplot.xlim([-0.1,5.3])
    matplotlib.pyplot.ylim([-15.2,8.4])

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    # 2.2. build distribution of transcript sets (all, no DETs, DET+ and DET-) per log10 expression units
    if timepoint != 'tp.1':
        transcriptSetsDistributionAnalyzer()

    print('\t completed analysis per time point.\n')

print('building figure for all time points...')
matplotlib.pyplot.plot(totalSetx,totalSety,'o',alpha=0.0333,mew=0,color='black')
matplotlib.pyplot.plot(totalHollowx,totalHollowy,'o',alpha=0.0333,mew=0,color='tan')
for i in range(len(regressions)):
    m=regressions[i][0]; c=regressions[i][1]
    floor=numpy.arange(0,highestExpression+0.1,0.1)
    e=list(m*numpy.array(floor)+c)
    matplotlib.pyplot.plot(floor,e,'-',lw=2,color=theColor[timepoints[i]])

matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.generalTrend.all.pdf')
matplotlib.pyplot.clf()

# 2.2. plot pattern from model
print('working on model building...')

figureName='figures/TE.model.pdf'
theColors=['red','orange','green','blue']

for i in range(len(sat)):
    px=sat[i][0]
    py=sat[i][1]
    matplotlib.pyplot.plot(px,py,'o',color=theColors[i])

matplotlib.pyplot.xlabel('mRNA [TPM]')
matplotlib.pyplot.ylabel('predicted footprint [TPM]')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()

# 3. build figure of TE (y-axis) and transcript length (x-axis)
print('building relationship between transcript length and TE...')

x=[]; y=[]
shortTranscript=10000; longTranscript=0

for element in lengthControl:
    x.append(element[0])
    y.append(element[1])
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='black')
if min(x) < shortTranscript:
    shortTranscript=min(x)
if max(x) > longTranscript:
    longTranscript=max(x)

x=[]; y=[]
for element in lengthControlHollow:
    x.append(element[0])
    y.append(element[1])
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='tan')
 
for i in range(len(lengthControlRegressions)):
    m=lengthControlRegressions[i][0]; c=lengthControlRegressions[i][1]
    floor=numpy.arange(shortTranscript,longTranscript+0.1,0.1)
    e=list(m*numpy.array(floor)+c)
    matplotlib.pyplot.plot(floor,e,'-',lw=2,color=theColor[timepoints[i]])

matplotlib.pyplot.xlabel('transcript length (log$_{10}$ bp)')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.length.control.pdf')
matplotlib.pyplot.clf()

# 4. build figure of TE (y-axis) and transcript half-life (x-axis)
print('building relationship between half-life...')

# 4.1. main analysis
x=[]; y=[]
fast=10000; slow=0

for element in degControl:
    x.append(element[0])
    y.append(element[1])
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='black')

if min(x) < fast:
    fast=min(x)
if max(x) > slow:
    slow=max(x)

x=[]; y=[]
for element in degControlHollow:
    x.append(element[0])
    y.append(element[1])
matplotlib.pyplot.plot(x,y,'o',alpha=0.0333,mew=0,color='tan')

for i in range(len(degControlRegressions)):
    m=degControlRegressions[i][0]; c=degControlRegressions[i][1]
    floor=numpy.arange(fast,slow+0.1,0.1)
    e=list(m*numpy.array(floor)+c)
    matplotlib.pyplot.plot(floor,e,'-',lw=2,color=theColor[timepoints[i]])

matplotlib.pyplot.xlim([0,25])

matplotlib.pyplot.xlabel('half-life (min)')
matplotlib.pyplot.ylabel('log$_{2}$ TE')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/TE.half-life.control.pdf')
matplotlib.pyplot.clf()

# 4.2. half-life/expression relationship
time=[]; exp=[]
low=10000; high=0

for element in degExpControl:
    if element[0] < 25:
        time.append(element[0])
        exp.append(element[1])
matplotlib.pyplot.plot(exp,time,'o',alpha=0.0333,mew=0,color='black')

if min(exp) < low:
    low=min(exp)
if max(exp) > high:
    high=max(exp)

for i in range(len(degExpRegressions)):
    m=degExpRegressions[i][0]; c=degExpRegressions[i][1]
    floor=numpy.arange(low,high+0.1,0.1)
    e=list(m*numpy.array(floor)+c)
    #matplotlib.pyplot.plot(floor,e,'-',lw=2,color=theColor[timepoints[i]])

matplotlib.pyplot.ylabel('half-life (min)')
matplotlib.pyplot.xlabel('log$_{10}$ TPM')
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figures/exp.half-life.control.pdf')
matplotlib.pyplot.clf()

# 5. final message
print('... all done. Bless!')
