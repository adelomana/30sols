###
### this script creates a heatmap of ribosomal protein expression
###

import sys,numpy,copy
import matplotlib,matplotlib.pyplot
import scipy,scipy.stats
import statsmodels,statsmodels.nonparametric,statsmodels.nonparametric.smoothers_lowess

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    this function creates a dictionary for expression values as
    expression[trna/rbf][ribo-pt gene name][timepoint][replicate]=value
    '''

    expression={}
    
    geneNames=[]
    timepoints=[]
    replicates=[]

    for sampleType in sampleTypes:
        expressionDataFile=expressionDataDir+'normalizedCounts.{}.csv'.format(sampleType)

        with open(expressionDataFile,'r') as f:

            firstLine=f.readline()
            header=firstLine.split(',')
            sampleNames=header[1:]
            sampleNames[-1]=sampleNames[-1].replace('\n','')

            for line in f:
                vector=line.split(',')

                # geneName
                geneName=vector[0]
                if geneName in riboPtNames:

                    # gene Names
                    if geneName not in geneNames:
                        geneNames.append(geneName)

                    for i in range(len(vector)-1):
                    
                        # timepoint
                        timepoint='tp.{}'.format(int(sampleNames[i].split('.')[-1]))
                        if timepoint not in timepoints:
                            timepoints.append(timepoint)

                        # replicate
                        replicate='rep.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
                        if replicate not in replicates:
                            replicates.append(replicate)

                        # value
                        value=float(vector[i+1])

                        # make sure keys exist
                        if sampleType not in expression.keys():
                            expression[sampleType]={}
                        if geneName not in expression[sampleType].keys():
                            expression[sampleType][geneName]={}
                        if timepoint not in expression[sampleType][geneName].keys():
                            expression[sampleType][geneName][timepoint]={}

                        expression[sampleType][geneName][timepoint][replicate]=value

    # sort variables
    sampleTypes.sort()
    geneNames.sort()
    timepoints.sort()
    replicates.sort()

    return expression,sampleTypes,timepoints,replicates

def riboPtNamesReader():

    '''
    this function reads the ribosomal protein names
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
expressionDataDir='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
scratchDir='/Volumes/omics4tb/alomana/scratch/'
theColors=['red','orange','green','blue']
sampleTypes=['trna','rbf']

# 1. read data
riboPtNames=riboPtNamesReader()
expression,sampleTypes,timepoints,replicates=expressionReader()

# 2. process data
# 2.1. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.2. plotting figures
allx=[]; ally=[]
for timepoint in timepoints:
    x=[]; y=[]
    for name in riboPtNames:
        valuesRNA=[]
        valuesRibo=[]
        for replicate in replicates:
            
            value=expression['trna'][name][timepoint][replicate]
            valuesRNA.append(value)

            value=expression['rbf'][name][timepoint][replicate]
            valuesRibo.append(value)

        averageRNA=numpy.mean(valuesRNA)
        averageRibo=numpy.mean(valuesRibo)

        x.append(averageRNA)
        y.append(averageRibo)

    # add to all time points variable
    for element in x:
        allx.append(element)
    for element in y:
        ally.append(element)
    print(len(allx),len(ally))

    # define group A and group B transcripts

    
    
    
    # 2.3. plotting dots
    theColor=theColors[int(timepoint[-1])-1]
    matplotlib.pyplot.plot(x,y,'o',alpha=0.5,mew=0,ms=8,color=theColor,label='TP {}'.format(timepoint[-1]))

# 2.4. compute linear regression
slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(allx,ally)
print('linear regression')
print('slope',slope)
print('intercept',intercept)
print('r_value',r_value)
print('pvalue',p_value)
print('std_err',std_err)
print()

# 2.5. compute and plot linear regression line
resolution=0.1
newx=numpy.arange(min(allx),max(allx),resolution)
newy=slope*newx+intercept
idx=numpy.where(newy>0)
matplotlib.pyplot.plot(newx[idx],newy[idx],lw=4,color='black')

description='R$^2$={:.2f}\np={:.2e}\na={:.2f}'.format(r_value**2,p_value,slope)
matplotlib.pyplot.text(7,9,description)

matplotlib.pyplot.xlabel('RNA-seq, log$_2$(normalized counts)')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$(normalized counts)')

matplotlib.pyplot.yticks([4,6,8,10])

#matplotlib.pyplot.xlim([1,18])
#matplotlib.pyplot.ylim([1,18])

matplotlib.pyplot.legend(markerscale=1.5)
        
figureName='figure.single.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()    
