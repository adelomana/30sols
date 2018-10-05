###
### this script performs analysis as Fig 1b of Schafer et al. (PMID, 26007203) using the data processed by kallisto/sleuth pipeline
###

import sys,os,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def colorAssigner(geneName,fcx,fcy):

    '''
    Assigning colors depending of significance and fold-change.
    '''

    # f1. coding significances
    significances=[]

    if changes[geneName][0][1] < 0.05:
        significances.append(1)
    else:
        significances.append(0)

    if changes[geneName][1][1] < 0.05:
        significances.append(1)
    else:
        significances.append(0)

    # f.2. selecting colors
    if significances == [1,1]:
        if fcx > 1 and fcy > 1:
            theColor='black'
        elif fcx < -1 and fcy < -1:
            theColor='black'
        elif fcx < -1 and fcy > 1:
            theColor='magenta'
        elif fcx > 1 and fcy < -1:
            theColor='cyan'
            print('observed cyan, exiting because not prepare to encounter this observation.')
            sys.exit()
        else:
            theColor='dubious'

    # orange and green
    elif significances == [0,1]:
        if  fcy > 1 and abs(fcx) < 1:
            theColor='orange'
        elif fcy < -1 and abs(fcx) < 1:
            theColor='green'
        else:
            theColor='dubious'

    # red and blue
    elif significances == [1,0]:
        if fcx < -1 and abs(fcy) < 1:
            theColor='red'
        elif fcx > 1 and abs(fcy) < 1:
            theColor='blue'
        else:
            theColor='dubious'
            
    # yellow
    elif significances == [0,0]:
        if abs(fcx) < 1 and abs(fcy) < 1:
            theColor='yellow'
        else:
            theColor='dubious'
    

    return theColor

def grapher(changes,flag):

    '''
    This function builds the figures for all genes that are not tRNAs
    '''

    # defining groups: 0, non-significant; 1, RBF significant; 2, RNA significant; 3, RBF and RNA significants
    positions={}
    positions['yellow']=[[],[]]
    positions['blue']=[[],[]]
    positions['red']=[[],[]]
    positions['orange']=[[],[]]
    positions['green']=[[],[]]
    positions['black']=[[],[]]

    yellowFile=open('results/results.yellow.txt','w')
    blueFile=open('results/results.blue.txt','w')
    redFile=open('results/results.red.txt','w')
    orangeFile=open('results/results.orange.txt','w')
    greenFile=open('results/results.green.txt','w')
    blackPlusFile=open('results/results.black.plus.txt','w')
    blackMinusFile=open('results/results.black.minus.txt','w')
    dubiousFile=open('results/results.dubious.txt','w')

    # selecting probing genes
    probingGenes=[]
    forbiddenGenes=['gene-VNGRS06350'] # new tRNA gene that does not have old annotation VNGt***
    
    for geneName in changes:
        
        oldName='-'
        try:
            oldName=synonyms[geneName]
        except:
            pass
        
        # selection
        if 'VNGt' not in oldName and geneName not in forbiddenGenes:
            probingGenes.append(geneName)
        else:
            print('bypassing {}/{}...'.format(geneName,oldName))
        
    for geneName in probingGenes:
        
        # compute fold-changes
        fcx=changes[geneName][0][0]
        fcy=changes[geneName][1][0]

        # assigning color
        theColor=colorAssigner(geneName,fcx,fcy)

        # filling up the dictionaries.
        if theColor != 'dubious':
            positions[theColor][0].append(fcx); positions[theColor][1].append(fcy)

        # associate to a synonym
        oldName='-'
        try:
            oldName=synonyms[geneName]
        except:
            pass

        # associate to a protein function
        proteinFunction='-'
        try:
            proteinFunction=proteinFunctions[geneName]
        except:
            pass

        # filling files
        if theColor == 'yellow':
            yellowFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'blue':
            blueFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'red':
            redFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'orange':
            orangeFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'green':
            greenFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'black':
            if fcx > 1:
                blackPlusFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
            else:
                blackMinusFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        elif theColor == 'dubious':
            dubiousFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(theColor,oldName,geneName,proteinFunction,fcx,fcy,changes[geneName][0][1],changes[geneName][1][1]))
        else:
            print('color not recognized')
            sys.exit()

    yellowFile.close()
    blueFile.close()
    redFile.close()
    orangeFile.close()
    greenFile.close()
    blackPlusFile.close()
    blackMinusFile.close()
    dubiousFile.close()
                
    # make figure for ribo-pt genes
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(1,1,1)
    majorTicks=numpy.arange(-10,10,5)
    minorTicks=numpy.arange(-10,10,1)

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
    matplotlib.pyplot.plot(positions['yellow'][0],positions['yellow'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='yellow',label='NR')
    matplotlib.pyplot.plot(positions['red'][0],positions['red'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='red',label='COMP+')
    matplotlib.pyplot.plot(positions['blue'][0],positions['blue'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='blue',label='COMP-')
    matplotlib.pyplot.plot(positions['orange'][0],positions['orange'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='orange',label='TL+')
    matplotlib.pyplot.plot(positions['green'][0],positions['green'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='green',label='TL-')
    matplotlib.pyplot.plot(positions['black'][0],positions['black'][1],'o',alpha=theAlpha,mew=0,ms=theSize,color='black',label='TR')

    matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
    matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

    matplotlib.pyplot.xlim([-10,10])
    matplotlib.pyplot.ylim([-10,10])

    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,ncol=2,fontsize=11)

    ax.grid(which='minor',alpha=0.2, ls=':')
    ax.grid(which='major',alpha=0.5,ls=':')

    ax.tick_params(which='minor',left='off',bottom='off') # turn off minor ticks

    matplotlib.pyplot.plot([-10,10],[-10,10],':',color='black')

    figureName='figure.{}.all.pdf'.format(flag)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

def significanceReaderHD():

    '''
    This function retrieves the significance defined by DESeq2.
    '''

    rawPositions={} # rawPositions[geneName]=[[trna.fc,trna.p-value],[rbf.fc,rbf.p-value]]

    sampleTypes=['trna','rbf']

    for sampleType in sampleTypes:
        dataFile=dataDirHD+'significance.{}.condition_tp.4_vs_tp.1.csv'.format(sampleType)
        
        with open(dataFile,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                
                geneName=(vector[0].replace('_','')).replace('"','')
                
                if vector[-1] != 'NA\n':
                    adj=float(vector[-1])
                else:
                    adj=1

                log2fc=float(vector[2])

                if geneName not in rawPositions:
                    rawPositions[geneName]=[[log2fc,adj]]
                else:
                    rawPositions[geneName].append([log2fc,adj])

    limited=[]
    for geneName in rawPositions:
        if len(rawPositions[geneName]) != 2:
            limited.append(geneName)
    changes={}
    for geneName in rawPositions:
        if geneName not in limited:
            changes[geneName]=rawPositions[geneName]
            
    return changes

def synonymsReader():

    '''
    This function reads a GFF3 file and returns synonyms. EG: VNG_RS00010 --> VNG0001H.
    '''

    synonyms={}
    definedGenes=[]
    proteinFunctions={}

    with open(gff3File,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1].replace('\n','')

                if len(vector) > 2:
                    if vector[2] == 'gene':
                        geneID=info.split('ID=')[1].split(';')[0]
                        definedGenes.append(geneID)
                
                if 'old_locus_tag=' in info:
                    old=info.split('old_locus_tag=')[1].split(';')[0]
                    new=info.split('ID=')[1].split(';')[0].replace('_','')

                    if '%' in old:
                        olds=old.split('%2C')
                        for element in olds:
                            synonyms[new]=element
                    else:
                        synonyms[new]=old

                # obtaining info about protein names
                if len(vector) > 2:
                    if vector[1] == 'Protein Homology':
                        geneID=info.split('Parent=')[1].split(';')[0].replace('_','')
                        proteinFunction=info.split('product=')[1].split(';')[0]
                        proteinFunctions[geneID]=proteinFunction

    definedGenes.sort()

    return synonyms,proteinFunctions

###
### MAIN
###

# 0. user-defined variables
gff3File='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
flag='hd'
dataDirHD='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'

# 1. obtain synonym
print('reading annotation file...')
synonyms,proteinFunctions=synonymsReader()

# 2. process data for STAR/htseq-count/DESeq2
print('reading RNA-seq and Ribo-seq info...')
changes=significanceReaderHD()

# 3. build figure
grapher(changes,flag)
