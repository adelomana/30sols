import sys
import matplotlib,matplotlib.pyplot
import numpy,numpy.linalg

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


# 0. user defined variables
theColor={}
theColor['tp.1']='red'
theColor['tp.2']='orange'
theColor['tp.3']='green'
theColor['tp.4']='blue'

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 2. computing the figure
print('computing the analysis...')

# 2.1. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. remove outliers
sat=[]

for timepoint in timepoints:

    figureName='figures/outfinder.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    
    for geneName in geneNames:
        
        # check consistency of mRNA
        mRNA_TPMs=[]; footprint_TPMs=[]
        for replicate in replicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate][timepoint][geneName])
            footprint_TPMs.append(rnaExpression['rbf'][replicate][timepoint][geneName])
            
        # data transformations 
        m=numpy.median(mRNA_TPMs); f=numpy.median(footprint_TPMs)
        
        log2m=numpy.log2(m+1)
        log10m=numpy.log10(m+1)
        log2f=numpy.log2(f+1)
        r=log2f-log2m

        if f == 0:
            hollowx.append(log10m); hollowy.append(r)
        else:
            setx.append(log10m); sety.append(r)

    # performing regression analysis
    A=numpy.vstack([setx,numpy.ones(len(setx))]).T
    solution,residuals,rank,s=numpy.linalg.lstsq(A,sety)

    m=solution[0]
    c=solution[1]
    expected=list(m*numpy.array(setx)+c)

    print(m,c)

    # computed from Matt Wall on log2
    ### satx=2**(numpy.array(setx))-1
    ### saty=(2**c)*((satx+1)**(1+m))-1

    # computed ALO on log10 for x
    satx=(10**numpy.array(setx))-1
    factor=m*numpy.log10(satx+1) + c + numpy.log2(satx+1)
    saty=(2**(factor))-1

    sat.append([list(satx),list(saty)])
        
    # plotting
    matplotlib.pyplot.plot(setx,sety,'o',alpha=0.1,mew=0,color='black')
    matplotlib.pyplot.plot(hollowx,hollowy,'o',alpha=0.1,mew=0,color='tab:brown')
    
    matplotlib.pyplot.plot(setx,expected,'-',lw=2,color=theColor[timepoint])

    print(numpy.max(setx))
    
    matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
    matplotlib.pyplot.ylabel('footprint/mRNA [log$_{2}$ ratio]')

    matplotlib.pyplot.xlim([-0.1,5.2])
    matplotlib.pyplot.ylim([-15,5.25])

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

# plotting the saturation effect
figureName='figures/saturation.pdf'
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
