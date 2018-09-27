###
### This script generates a scatter plot for ribosome loading across expression
### 

import sys,math,pandas,seaborn
import matplotlib,matplotlib.pyplot
import numpy,numpy.linalg
import scipy,scipy.stats

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

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

# 2. perform analysis
print('computing the analysis...')

# 2.0. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. build figures on general pattern
sat=[] # needed for second part of the script, the pattern model

totalSetx=[]; totalSety=[]; totalHollowx=[]; totalHollowy=[]

for timepoint in timepoints:

    figureName='figures/panel.A.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    cloudColors=[]; groundColors=[]
    
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
                
        else:
            if rsem_mRNA < 0.3 and rsem_RF < 0.3:        

                setx.append(m); sety.append(r)
                totalSetx.append(m); totalSety.append(r)                

               
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

    # computed from Matt Wall on log2
    ### satx=2**(numpy.array(setx))-1
    ### saty=(2**c)*((satx+1)**(1+m))-1

    # computed ALO on log10 for x #!!! check model
    #satx=(10**numpy.array(setx))-1
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
matplotlib.pyplot.savefig('figures/TE.generalTrend.all.pdf')
matplotlib.pyplot.clf()

# 3. plot pattern from model
print('working on model building...')

figureName='figures/TE.model.pdf'
theColors=['red','orange','green','blue']

for i in range(len(sat)):
    px=sat[i][0]
    py=sat[i][1]
    matplotlib.pyplot.plot(px,py,'-',lw=2,color=theColors[i])
    print('last value',py[-1])

matplotlib.pyplot.xlabel('mRNA [TPM]')
matplotlib.pyplot.ylabel('predicted footprint [TPM]')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()


