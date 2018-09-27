###
### This script builds a scatter plot of mRNA vs RF (one per order of magnitude of mRNA) for up regulated proteins.
###

import sys,os,numpy
import scipy,scipy.stats

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def DETreader():

    '''
    This function reads the DETs from sleuth and provides a dictionary as 
    DETs[timepoint][geneName]=1/-1
    DETs of abs(FC) < 2 and max expression < 10 TPMs will be excluded.
    '''

    DETs={}
    # iterate over each time point to retrieve the DETs
    for timepoint in rnaTimepoints[1:]:
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
                rna0=numpy.mean([rnaExpression['trna'][replicate]['tp.1'][geneName] for replicate in rnaReplicates])
                rna1=numpy.mean([rnaExpression['trna'][replicate][timepoint][geneName] for replicate in rnaReplicates])
                log2fc=numpy.log2(rna1/rna0)

                if abs(log2fc) > 1 and numpy.max([rna0,rna1]) > 10:
                    if log2fc > 0:
                        flag=1
                    else:
                        flag=-1
                    DETs[timepoint][geneName]=flag

    return DETs

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

def proteomicsDataReader():

    '''
    This function reads data available and outputs the defined dictionary.
    '''

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=[log2FC,pvalue]

    conditions=[]
    geneNames=[]
    timepoints=[]
    replicates=[]
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    lostGenes=[]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in conditions:
            conditions.append(condition)
        if replicate not in replicates:
            replicates.append(replicate)

        if condition not in data.keys():
            data[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}

        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')

                try:
                    geneName=synonyms[vector[0]]
                    geneName=geneName.replace('gene-','').replace('_','')

                    if geneName not in geneNames:
                        geneNames.append(geneName)

                    a=float(vector[2])
                    b=float(vector[6])
                    c=float(vector[10])

                    d=float(vector[4])
                    e=float(vector[8])
                    f=float(vector[12])

                    data[condition][replicate]['tp2vs1'][geneName]=[a,d]
                    data[condition][replicate]['tp3vs1'][geneName]=[b,e]
                    data[condition][replicate]['tp4vs1'][geneName]=[c,f]
                    
                except:
                    if vector[0] not in lostGenes:
                        lostGenes.append(vector[0])
                
    # sort
    geneNames.sort()
    conditions.sort()
    replicates.sort()

    # print lost genes
    print('\t {} protein lost because of new annotation.'.format(len(lostGenes)))
    
    return data,geneNames,conditions,replicates,timepoints

def scatterPlotBuilder():

    '''
    This function builds a mRNA/footprints scatterplot based on FC.
    '''

    return None

def synonymsReader():

    '''
    This function reads the GFF3 file and returns a dictionary with synonyms between old and new locus names.
    '''

    synonyms={}
    with open(annotationFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1].replace('\n','')
                if 'old_locus_tag=' in info:
                    old=info.split('old_locus_tag=')[1].split(';')[0]
                    new=info.split('ID=')[1].split(';')[0]

                    if '%' in old:
                        olds=old.split('%2C')
                        for element in olds:
                            synonyms[element]=new
                    else:
                        synonyms[old]=new

    # reverting mapping
    synonymsReverseMapping={v: k for k, v in synonyms.items()}
                
    return synonyms,synonymsReverseMapping

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
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
transcriptomeAnnotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.1.cs.NC_002608.1.fasta'
transcriptomicsDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'
DETsDir='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/'

alpha=0.05

slope=-1.26587161237
intercept=3.32469140168

# 1. reading data
print('reading data...')

# 1.1. define synonyms
synonyms,synonymsReverseMapping=synonymsReader()
NCsynonyms=NCsynonymsReader()

# 1.1. reading mRNA data
rnaExpression,rnaNames,rnaTimepoints,rnaReplicates=transcriptomicsReader()
DETs=DETreader()

# 1.2. read proteomics data
print('\t reading proteomics data...')
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteomicsDataReader()

# 2. analysis
print('performing analysis...')

# 2.1. get up-regulated proteins
upPts=[]
for condition in proteinConditions:
    figureFile='figure.up.{}.pdf'.format(condition)
    print(condition)

    # 2.1. get up-regulated proteins
    for proteinName in proteinNames:
        try:
            pvalues=[proteinAbundance[condition][replicate]['tp4vs1'][proteinName][1] for replicate in proteinReplicates]
        except:
            pass
            #print('\t loosing {} because failed to be detected at this time point.'.format(proteinName))
        if numpy.max(pvalues) < alpha:
            upPts.append(proteinName)
    print('\t {} proteins found upregulated.'.format(len(upPts)))

    # 2.2. filter out genes that are transcriptionally upregulated
    candidates=[]; 
    for protein in upPts:
        if protein not in list(DETs['tp.4'].keys()):
            candidates.append(protein)
    print('\t {} upregulated proteins do not have mRNA up-regulated.'.format(len(candidates)))

    # 2.3. filter out transcripts that are highly expressed: > 100 TPMs or they have noise
    finalSet=[]
    for name in candidates:

        # check consistency of mRNA
        mRNA_TPMs=[]
        for replicate in rnaReplicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate]['tp.4'][name])
            
        # data transformations and quality check
        log2M=numpy.log2(numpy.array(mRNA_TPMs)+1)

        # noise
        if numpy.max(log2M) > numpy.log2(11): # if expression is below 10 TPMs, don't consider noise
            sem=numpy.std(log2M)/numpy.sqrt(len(log2M))
            rsem_mRNA=sem/numpy.mean(log2M)
        else:
            rsem_RNA=0

        # condition
        if rsem_mRNA < 0.3 and numpy.median(log2M) < 5:
            finalSet.append(name)
    print('\t {} low abundance transcripts.'.format(len(finalSet)))            
    
    # 2.4. build a scatterplot
    setx=[]; sety=[]
    hollowx=[]; hollowy=[]
    for name in finalSet:

        # x value face value of mRNA
        mRNA_TPMs=[]; footprint_TPMs=[]
        for replicate in rnaReplicates:
            mRNA_TPMs.append(rnaExpression['trna'][replicate]['tp.4'][name])
            footprint_TPMs.append(rnaExpression['rbf'][replicate]['tp.4'][name])
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
        else:
            if rsem_mRNA < 0.3 and rsem_RF < 0.3:        
                setx.append(m); sety.append(r)
        

    # plot figure
    for i in range(len(setx)):
        matplotlib.pyplot.plot(setx[i],sety[i],'o',alpha=0.5,mew=0,color='black')
    for i in range(len(hollowx)):
        matplotlib.pyplot.plot(hollowx[i],hollowy[i],'o',alpha=0.5,mew=0,color='black')

    # plot line
    m=slope
    c=intercept
    expected=list(m*numpy.array(setx)+c)
    matplotlib.pyplot.plot(setx,expected,'-',lw=1,color='black')

    #matplotlib.pyplot.xlim([-0.1,5.3])
    #matplotlib.pyplot.ylim([-15.2,8.4])
    
    matplotlib.pyplot.xlabel('mRNA [log$_{10}$ TPM+1]')
    matplotlib.pyplot.ylabel('footprint/mRNA [log$_{2}$ ratio]')

    matplotlib.pyplot.grid(True,alpha=0.5,ls=':')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()
    
    #
    print()

# for each order of magnitud of expression (median of t4 and t1), build a FC of mRNA vs RF in the context of the model line
