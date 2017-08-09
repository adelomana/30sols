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

# 0.1. paths
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression/expressionMatrix.kallisto.txt'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 2. computing the figure
print('computing the analysis...')

# 2.1. remove outliers

sat=[]

for timepoint in timepoints:

    figureName='figures/outfinder.{}.pdf'.format(timepoint)

    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    
    for geneName in geneNames:
        
        # check consistency of mRNA
        mRNA_TPM=[]; footprint_TPM=[]
        for replicate in replicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate][timepoint][geneName])
            footprint_TPMs.append(rnaExpression['rbf'][replicate][timepoint][geneName])
            
        # data transformations 
            
        e=numpy.log2(x+1)
        l=numpy.log10(x+1)
        f=numpy.log2(y+1)
        r=f-e

            if y == 0:
                hollowx.append(l); hollowy.append(r)
            else:
                setx.append(l); sety.append(r)

    # performing regression analysis
    A=numpy.vstack([setx,numpy.ones(len(setx))]).T
    solution,residuals,rank,s=numpy.linalg.lstsq(A,sety)

    m=solution[0]
    c=solution[1]
    expected=list(m*numpy.array(setx)+c)

    print(m,c)

    satx=10**(numpy.array(setx))-1
    saty=(10**c)*(satx**(1+m))-1
    sat.append([satx,saty])
        
    # plotting
    matplotlib.pyplot.plot(setx,sety,'o',alpha=0.1,mew=0,color='black')
    matplotlib.pyplot.plot(hollowx,hollowy,'o',alpha=0.1,mew=0,color='orange')
    
    matplotlib.pyplot.plot(setx,expected,'-',lw=2,color='green')
    
    matplotlib.pyplot.xlabel('mRNA [log$_2$ TPM+1]')
    matplotlib.pyplot.ylabel('footprint/mRNA [log$_{10}$ FC]')

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()


# plotting the saturation effect
figureName='figures/saturation.pdf'
theColors=['blue','green','orange','red']
for i in range(len(sat)):
    matplotlib.pyplot.semilogx(sat[i][0],sat[i][1],'-',lw=2,color=theColors[i])

matplotlib.pyplot.xlabel('mRNA [TPM]')
matplotlib.pyplot.ylabel('predicted footprint [TPM]')
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()
