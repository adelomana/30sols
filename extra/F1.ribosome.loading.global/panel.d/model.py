###
### this script performs analysis as Fig 1b of Schafer et al. (PMID, 26007203) using the data processed by kallisto/sleuth pipeline
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def colorAssigner(geneName,fcx,fcy):

    """
    code 0: no-significant abs(FC) < 1
    """

    code=0

    if geneName in significantRNAs:
        if fcx > 1:
            code=1
        elif fcx < -1:
            code=2
            
    if geneName in significantRBFs:
        if fcy > 1:
            code=code+10
        elif fcy < -1:
            code=code+20

    if code == 0:
        theColor='yellow'
    elif code == 1:
        theColor='blue'
    elif code == 2:
        theColor='red'
    elif code == 10:
        theColor='orange'
    elif code == 20:
        theColor='green'
    elif code == 21:
        theColor='magenta'
    elif code == 12:
        theColor='magenta'
    else:
        theColor='black'

    return theColor

def expressionReaderHD():

    """
    idem as expressionReaderKS but reading data from STAR/htseq-count/DESeq2 pipeline
    """

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    sampleTypes=['trna','rbf']

    for sampleType in sampleTypes:
        dataFile=dataDirHD+'normalizedCounts.{}.csv'.format(sampleType)

        with open(dataFile,'r') as f:

            firstLine=f.readline()
            header=firstLine.split(',')
            sampleNames=header[1:]
            sampleNames[-1]=(sampleNames[-1].replace('\n','')).replace('"','')
            sampleNames=[element.replace('"','') for element in sampleNames]

            for line in f:
                vector=line.split(',')

                # geneName
                geneName=(vector[0].replace('_','')).replace('"','')
                if geneName not in geneNames:
                    geneNames.append(geneName)
                    
                for i in range(len(sampleNames)):

                    # sampleType
                    if sampleType not in sampleTypes:
                        sampleTypes.append(sampleType)

                    # timepoint
                    timepoint='timepoint.{}'.format(int(sampleNames[i].split('.')[-2]))
                    if timepoint not in timepoints:
                        timepoints.append(timepoint)

                    # replicate
                    replicate='replicate.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
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

    return expression,sampleTypes,geneNames,timepoints,replicates

def expressionReaderSC():

    """
    idem as expressionReaderKS but reading data from STAR/cuffdiff pipeline
    """

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    with open(expressionDataFileSC,'r') as f:

        firstLine=f.readline()
        header=firstLine.split('\t')
        sampleNames=header[1:]
        sampleNames[-1]=sampleNames[-1].replace('\n','')
        sampleNames=[element.replace('_0','') for element in sampleNames]

        for line in f:
            vector=line.split('\t')

            # geneName
            geneName=vector[0].replace('_','')            
            if geneName not in geneNames:
                geneNames.append(geneName)

            for i in range(len(sampleNames)):

                # sampleType
                sampleType=sampleNames[i].split('.')[0]
                if sampleType not in sampleTypes:
                    sampleTypes.append(sampleType)

                # timepoint
                timepoint='timepoint.{}'.format(int(sampleNames[i][-1]))
                if timepoint not in timepoints:
                    timepoints.append(timepoint)

                # replicate
                replicate='replicate.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
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
    
    return expression,sampleTypes,geneNames,timepoints,replicates

def expressionReaderKS():

    '''
    this function creates a dictionary for expression values as
    expression[trna/rbf][geneName][timepoint][replicate]=value
    '''

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    with open(expressionDataFileKS,'r') as f:

        firstLine=f.readline()
        header=firstLine.split('\t')
        sampleNames=header[1:-1]

        for line in f:
            vector=line.split('\t')

            # geneName
            geneName=vector[0].replace('_','')            
            if geneName not in geneNames:
                geneNames.append(geneName)

            for i in range(len(sampleNames)):

                # sampleType
                sampleType=sampleNames[i].split('.')[0]
                if sampleType not in sampleTypes:
                    sampleTypes.append(sampleType)

                # timepoint
                timepoint='timepoint.{}'.format(int(sampleNames[i][-1]))
                if timepoint not in timepoints:
                    timepoints.append(timepoint)

                # replicate
                replicate='replicate.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
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
    
    return expression,sampleTypes,geneNames,timepoints,replicates

def grapher(flag):

    '''
    this function builds the figures for all or ribosomal protein genes
    '''

    # defining groups: 0, non-significant; 1, RBF significant; 2, RNA significant; 3, RBF and RNA significants
    positions={}
    positions['ribo.yellow']=[[],[]]; positions['all.yellow']=[[],[]]
    positions['ribo.blue']=[[],[]]; positions['all.blue']=[[],[]]
    positions['ribo.red']=[[],[]]; positions['all.red']=[[],[]]
    positions['ribo.orange']=[[],[]]; positions['all.orange']=[[],[]]
    positions['ribo.green']=[[],[]]; positions['all.green']=[[],[]]
    positions['ribo.magenta']=[[],[]]; positions['all.magenta']=[[],[]]
    positions['ribo.black']=[[],[]]; positions['all.black']=[[],[]]
   
    
    timepointLate=timepoints[-1]
    timepointEarly=timepoints[0]
     
    for geneName in geneNames:

        # compute averages for RNA-seq late
        x=numpy.mean([numpy.log2(expression[sampleTypes[1]][geneName][timepointLate][replicate]+1) for replicate in replicates])

        # compute averages for RNA-seq early
        y=numpy.mean([numpy.log2(expression[sampleTypes[1]][geneName][timepointEarly][replicate]+1) for replicate in replicates])

        try: # counts in Ribo-seq below sum(10) are removed
            # compute averages for Ribo-seq late
            z=numpy.mean([numpy.log2(expression[sampleTypes[0]][geneName][timepointLate][replicate]+1) for replicate in replicates])

            # compute averages for Ribo-seq early
            w=numpy.mean([numpy.log2(expression[sampleTypes[0]][geneName][timepointEarly][replicate]+1) for replicate in replicates])
        except:
            z=0; w=0
            if geneName in riboPtNames:
                extra=1
            else:
                extra=0
            #print('WARNING: {} ({}) has very low counts (sum <10) for Ribo-seq experiment.'.format(geneName,extra))
        
        # compute fold-changes
        fcx=x-y
        fcy=z-w

        # assigning color
        theColor=colorAssigner(geneName,fcx,fcy)

        # filling up the dictionaries.
        tag='all.'+theColor
        positions[tag][0].append(fcx); positions[tag][1].append(fcy)

        if geneName in riboPtNames:
            tag='ribo.'+theColor
            positions[tag][0].append(fcx); positions[tag][1].append(fcy)
                
    # make figure for ribo-pt genes
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(1,1,1)
    majorTicks=numpy.arange(-10,10,5)
    minorTicks=numpy.arange(-10,10,1)

    ax.set_xticks(majorTicks)
    ax.set_xticks(minorTicks,minor=True)
    ax.set_yticks(majorTicks)
    ax.set_yticks(minorTicks,minor=True)

    theAlpha=0.3
    theSize=8
    matplotlib.pyplot.plot(positions['ribo.yellow'][0],positions['ribo.yellow'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='yellow',label='NR')
    matplotlib.pyplot.plot(positions['ribo.red'][0],positions['ribo.red'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='red',label='COMP+')
    matplotlib.pyplot.plot(positions['ribo.blue'][0],positions['ribo.blue'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='blue',label='COMP-')
    matplotlib.pyplot.plot(positions['ribo.orange'][0],positions['ribo.orange'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='orange',label='TL+')
    matplotlib.pyplot.plot(positions['ribo.green'][0],positions['ribo.green'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='green',label='TL-')
    matplotlib.pyplot.plot(positions['ribo.magenta'][0],positions['ribo.magenta'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='magenta',label='OPP')
    matplotlib.pyplot.plot(positions['ribo.black'][0],positions['ribo.black'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='black',label='TR')
   
    matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
    matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

    matplotlib.pyplot.xlim([-10,10])
    matplotlib.pyplot.ylim([-10,10])

    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,ncol=2,fontsize=11)

    ax.grid(which='minor',alpha=0.2, ls=':')
    ax.grid(which='major',alpha=0.5,ls=':')

    ax.tick_params(which='minor',left='off',bottom='off') # turn off minor ticks

    matplotlib.pyplot.plot([-10,10],[-10,10],':',color='black')

    figureName='figure.{}.ribo.pdf'.format(flag)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    # make figure for all genes
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(1,1,1)
    majorTicks=numpy.arange(-10,10,5)
    minorTicks=numpy.arange(-10,10,1)

    ax.set_xticks(majorTicks)
    ax.set_xticks(minorTicks,minor=True)
    ax.set_yticks(majorTicks)
    ax.set_yticks(minorTicks,minor=True)

    theAlpha=0.3
    theSize=8
    matplotlib.pyplot.plot(positions['all.yellow'][0],positions['all.yellow'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='yellow',label='NR')
    matplotlib.pyplot.plot(positions['all.red'][0],positions['all.red'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='red',label='COMP+')
    matplotlib.pyplot.plot(positions['all.blue'][0],positions['all.blue'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='blue',label='COMP-')
    matplotlib.pyplot.plot(positions['all.orange'][0],positions['all.orange'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='orange',label='TL+')
    matplotlib.pyplot.plot(positions['all.green'][0],positions['all.green'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='green',label='TL-')
    matplotlib.pyplot.plot(positions['all.magenta'][0],positions['all.magenta'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='magenta',label='OPP')
    matplotlib.pyplot.plot(positions['all.black'][0],positions['all.black'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='black',label='TR')

    matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
    matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

    matplotlib.pyplot.xlim([-10,10])
    matplotlib.pyplot.ylim([-10,10])

    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,ncol=2,fontsize=11)

    ax.grid(which='minor',alpha=0.2, ls=':')
    ax.grid(which='major',alpha=0.5,ls=':')

    ax.tick_params(which='minor',left='off',bottom='off') # turn off minor ticks

    matplotlib.pyplot.plot([-10,10],[-10,10],':',color='black')

    figureName='figure.{}.all.png'.format(flag)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

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

def significanceReaderHD():

    '''
    this function retrieves the significance defined by DESeq2
    '''

    significantRBFs=[]; significantRNAs=[]

    sampleTypes=['trna','rbf']

    for sampleType in sampleTypes:
        dataFile=dataDirHD+'significance.{}.csv'.format(sampleType)
        
        with open(dataFile,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                geneName=(vector[0].replace('_','')).replace('"','')
                
                if vector[-1] != 'NA\n':
                    adj=float(vector[-1])
                else:
                    adj=1

                if adj < 0.05:
                    if sampleType == 'trna':
                        significantRNAs.append(geneName)
                    elif sampleType == 'rbf':
                        significantRBFs.append(geneName)
                    else:
                        print('error 361')
                        sys.exit()

    return significantRBFs,significantRNAs

def significanceReaderKS():

    '''
    This function retrieves significance defined by Sleuth.
    '''

    significantRBFs=[]; significantRNAs=[]

    workingFiles=os.listdir(significanceDirKS)
    
    for workingFile in workingFiles:
        with open(significanceDirKS+workingFile,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                target=vector[1].replace('"','')


                print(significanceDirKS+workingFile)
                print(target)
                
                geneName=synonyms[target].replace('_','')
                if 'RBF' in workingFile:
                    significantRBFs.append(geneName)
                elif 'RNA' in workingFile:
                    significantRNAs.append(geneName)
                else:
                    print('error a')
                    sys.exit()

    return significantRBFs,significantRNAs

def significanceReaderSC():

    '''
    this function retrieves the significance defined by Cuffdiff
    '''

    significantRBFs=[]; significantRNAs=[]

    allDirs=os.listdir(significanceDirSC)
    workingDirs=[element for element in allDirs if 'vs' in element]
    
    for workingDir in workingDirs:
        workingFile=significanceDirSC+workingDir+'/gene_exp.diff'
        with open(workingFile,'r') as f:
            next(f)
            for line in f:
                vector=line.split('\t')
                geneName=vector[0].replace('_','')
                significance=vector[-1].replace('\n','')
                if significance == 'yes':
                    if 'RBF' in workingFile:
                        significantRBFs.append(geneName)
                    elif 'RNA' in workingFile:
                        significantRNAs.append(geneName)
                    else:
                        print('error a')
                        sys.exit()

    return significantRBFs,significantRNAs

def synonymReader():

    synonyms={}
    with open(transcriptomeFastaFile,'r') as f:
        for line in f:
            if line[0] == '>':

                id=line.split('>')[1].split(' ')[0]
                if 'locus_tag' in line:
                    name=line.split('[locus_tag=')[1].split(' ')[0].replace(']','')
                elif 'gene' in line:
                    name=line.split('[gene=')[1].split(' ')[0].replace(']','')
                else:
                    print('error while parsing transcriptome fasta file...')
                    sys.exit()

                synonyms[id]=name

    return synonyms

### MAIN

# 1. read general data
print('reading data...')
transcriptomeFastaFile='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/hsa.ASM680v1.fasta'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
synonyms=synonymReader()
riboPtNames=riboPtNamesReader()

# 2. process data for kallisto/sleuth
print('')
print('processing kallisto/sleuth...')

flag='ks'

expressionDataFileKS='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'
significanceDirKS='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/'

expression,sampleTypes,geneNames,timepoints,replicates=expressionReaderKS()
significantRBFs,significantRNAs=significanceReaderKS()
print('genes retrieved {}'.format(len(geneNames)))
print('number of significant transcripts {}'.format(len(significantRNAs)))
print('number of significant footprints {}'.format(len(significantRBFs)))
grapher(flag)

# 3. process data for STAR/cuffdiff
print('')
print('processing STAR/cuffdiff...')

flag='sc'

expressionDataFileSC='/Volumes/omics4tb/alomana/projects/TLR/data/cufflinks/allSamples/genes.fpkm_table'
significanceDirSC='/Volumes/omics4tb/alomana/projects/TLR/data/cufflinks/'

expression,sampleTypes,geneNames,timepoints,replicates=expressionReaderSC()
significantRBFs,significantRNAs=significanceReaderSC()
print('genes retrieved {}'.format(len(geneNames)))
print('number of significant transcripts {}'.format(len(significantRNAs)))
print('number of significant footprints {}'.format(len(significantRBFs)))
grapher(flag)

# 3. process data for STAR/htseq-count/DESeq2
print('')
print('processing STAR/htseq-count/DESeq2...')

flag='hd'
dataDirHD='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'

expression,sampleTypes,geneNames,timepoints,replicates=expressionReaderHD()
significantRBFs,significantRNAs=significanceReaderHD()
print('genes retrieved {}'.format(len(geneNames)))
print('number of significant transcripts {}'.format(len(significantRNAs)))
print('number of significant footprints {}'.format(len(significantRBFs)))
grapher(flag)
