
"""
iterate over timepoints
iterate over genes
if both are below one, discard
sum +1, convert to log2 or log10
check that there is replicates heterogeneity
plot
"""

import os,sys,numpy
import scipy,scipy.stats
import matplotlib, matplotlib.pyplot

def dataReader():

    with open(expressionDataFile,'r') as f:

        # dealing with the header
        header=f.readline()
        vectorHeader=header.split('\t')[1:]
        experiments=[]; timepoints=[]; replicates=[]
        for element in vectorHeader:

            if element.split('.')[0] == 'trna':
                experiments.append('mRNA')
            elif element.split('.')[0] == 'rbf':
                experiments.append('footprint')
            else:
                print('error 1')
                sys.exit()

            timepoints.append(int(element.split('.')[-1].split('_')[0]))

            replicates.append(int(element.split('rep.')[1].split('.tp')[0]))

        # filling up the expression dictionary: expression[mRNA/footprint][tp#][rep#][geneName]=float
        expression={}
        for experiment in set(experiments):
            expression[experiment]={}
            for timepoint in set(timepoints):
                expression[experiment][timepoint]={}
                for replicate in replicates:
                    expression[experiment][timepoint][replicate]={}
                    
        # dealing with the body of expression
        geneNames=[]
        for line in f:
            vector=line.split('\t')
            geneName=vector[0]
            geneNames.append(geneName)
            measurements=[float(element) for element in vector[1:]]
            for i in range(len(measurements)):
                expression[experiments[i]][timepoints[i]][replicates[i]][geneName]=measurements[i]            
        
    return expression,experiments,timepoints,replicates,geneNames

# 0. user defined variables
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/cufflinks/allSamples/genes.fpkm_table'
noiseThreshold=2

# 1. reading the data
print('reading expression...')
expression,experiments,timepoints,replicates,geneNames=dataReader()
print(len(geneNames))

# 2. computing correlations
print('computing correlations...')


variationT=[]
variationF=[]
          
for timepoint in set(timepoints):

    figureFileName='figure.timepoint.{}.pdf'.format(timepoint)
    x=[]
    y=[]
    
    for geneName in geneNames:
        
        transcriptValues=[]; footprintValues=[]
        for replicate in set(replicates):
            transcriptValues.append(expression['mRNA'][timepoint][replicate][geneName])
            footprintValues.append(expression['footprint'][timepoint][replicate][geneName])

        # checking expression above noise
        if numpy.mean(transcriptValues) > 2:
            
            # checking for coherent expression
            logT=numpy.log10([value+1 for value in transcriptValues])
            logF=numpy.log10([value+1 for value in footprintValues])

            cvT=scipy.stats.variation(logT)
            if max(logF) > 0:
                cvF=scipy.stats.variation(logF)
            else:
                cvF=0

            #variationT.append(cvT)
            #variationF.append(cvF)

            if max([cvT,cvF]) < 0.2:

                x.append(numpy.mean(logT))
                y.append(numpy.mean(logF))

    print(len(x))
    matplotlib.pyplot.plot(x,y,'ok',alpha=0.2)
    matplotlib.pyplot.savefig(figureFileName)
    matplotlib.pyplot.clf()


#matplotlib.pyplot.hist(variationT)
#matplotlib.pyplot.savefig('figure.pdf')

            
# 2.1. cleaning low similarity replicates

# 2.2. plotting correlations
