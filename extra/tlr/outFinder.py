import sys
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

                # filter out abs(FC) < 2 or max expression < 10 TPMs
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

def NCsynonymsReader():

    '''
    This function builds a dictionary on the synonyms from NC to RS.
    '''

    NCsynonyms={}
    with open(transcriptomeAnnotationFile,'r') as f:
        for line in f:
            if line[0] == '>':
                vector=line.split(' ')
                ncName=vector[0].replace('>','')
                for element in vector:
                    if 'locus_tag' in element:
                        rsName=element.split('locus_tag=')[1].replace(']','').replace('_','')
                NCsynonyms[ncName]=rsName

    return NCsynonyms

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
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. reading data
print('reading data...')

# 1.1. reading mRNA data
rnaExpression,geneNames,timepoints,replicates=transcriptomicsReader()

# 1.2. read DETs
NCsynonyms=NCsynonymsReader()
DETs=DETreader()

# 2. perform analysis
print('computing the analysis...')

# 2.0. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# 2.1. build figures on general pattern
sat=[] # needed for second part of the script, the pattern model

for timepoint in timepoints:

    figureName='figures/outfinder.new.{}.pdf'.format(timepoint)

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
        if numpy.max(log2M) > numpy.log2(11): # if expression is below 10 TPMs, don't consider noise
            sem=numpy.std(log2M)/numpy.sqrt(len(log2M))
            rsem_mRNA=sem/numpy.mean(log2M)
        else:
            rsem_RNA=0
            
        if numpy.max(log2F) > numpy.log2(11): # if expression is below 10 TPMs, don't consider noise
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
                dotColor=dotColorFinder(timepoint,geneName,'ground')
                groundColors.append(dotColor)
        else:
            if rsem_mRNA < 0.3 and rsem_RF < 0.3:        
                setx.append(m); sety.append(r)
                dotColor=dotColorFinder(timepoint,geneName,'cloud')
                cloudColors.append(dotColor)

    # perform regression analysis
    print('\t regression results...')
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(setx,sety)
    print('\t linear regression')
    print('\t slope',slope)
    print('\t intercept',intercept)
    print('\t r_value',r_value)
    print('\t pvalue',p_value)
    print('\t std_err',std_err)
    print()

    # compute for the model
    m=slope
    c=intercept
    expected=list(m*numpy.array(setx)+c)

    # computed from Matt Wall on log2
    ### satx=2**(numpy.array(setx))-1
    ### saty=(2**c)*((satx+1)**(1+m))-1

    # computed ALO on log10 for x
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

    #matplotlib.pyplot.title('n = {}'.format(len(setx)))

    matplotlib.pyplot.xlim([-0.1,5.3])
    matplotlib.pyplot.ylim([-15.2,8.4])

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

# 2.2. plot pattern from model
figureName='figures/saturation.new.pdf'
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
