###
### This script generates a filtered set of DETs based on passing the rule abs log2 FC > 1
###

import os,pandas,sys,numpy

def formatter():

    '''
    This function generates a new file with only those DETs that pass the rule of abs log2 FC > 1
    '''

    label=sleuthResultsFile.replace('.csv','')

    # f.0. open comparisons files
    outputFile=sleuthResultsFile.replace('.csv','.filtered.txt')

    g=open(outputFile,'w')
    g.write('geneName\tmeanA\tmeanB\tlog2FC\tq-value\tq<0.05|abs(log2(FC))>1|maxExp>10.\n')

    formattedNum=[]; formattedDen=[]
    for sampleName in expression:
        if 'trna' in sampleName:
            if 'tp.4' in sampleName:
                formattedNum.append(sampleName)
            if 'tp.2' in sampleName:
                formattedDen.append(sampleName)

    # f.1. retain the number of passed
    passed={}   
    with open(sleuthResultsFile,'r') as f:
        next(f)
        for line in f:
            v=line.split(',')
            locusID=v[1].replace('"','')
            qValue=v[3]
            geneName=synonyms[locusID]

            # f.2. check for FC
            a=[expression[sample][geneName] for sample in formattedNum]
            b=[expression[sample][geneName] for sample in formattedDen]
            meanA=numpy.mean(a)
            meanB=numpy.mean(b)
            FC=(meanA+1)/(meanB+1)
            log2FC=numpy.log2(FC)

            if abs(log2FC) > 1 and numpy.max([meanA,meanB]) > 10.:
                flag=True
                passed[geneName]=log2FC
                g.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneName,meanA,meanB,log2FC,qValue,flag))
            else:
                flag=False
                print(geneName,log2FC,meanA,meanB,qValue)

            # f.3. write line in filtered file
            g.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneName,meanA,meanB,log2FC,qValue,flag))

    # f.4. close filtred file
    g.close()

    # f.5. last message
    print('\t {} DETs passed filter for comparison {}.'.format(len(passed),label))

    return passed

def expressionReader():

    '''
    This function creates a dictionary with the expression from kallisto.
    '''

    expression={}
    with open(expressionFile,'r') as f:
        
        h=f.readline()
        headers=h.split()
        for sample in headers:
            expression[sample]={}

        for line in f:
            v=line.split('\t')
            geneName=v[0]
            for i in range(len(v[2:])):
                value=float(v[i+1])
                expression[headers[i]][geneName]=value

    return expression

def synonymsFinder():

    '''
    This function builds a dictionary of gene-locus associations.
    '''
    
    synonyms={}
    with open(transcriptomeFastaFile,'r') as f:
        for line in f:
            if line[0] == '>':

                id=line.split('>')[1].split(' ')[0]
                if 'locus_tag' in line:
                    name=line.split('[locus_tag=')[1].split(' ')[0].replace(']','')
                else:
                    print('error while parsing transcriptome fasta file...')
                    sys.exit()

                synonyms[id]=name

    return synonyms

# 0. user defined variables
sleuthResultsFile='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/sleuthResultsRNA.41.csv'
expressionFile='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'
transcriptomeFastaFile='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.1.cs.NC_002608.1.fasta'

# 1. read data
print('reading data...')

# 1.3. define synonyms
synonyms=synonymsFinder()

# 1.4. define expression
expression=expressionReader()

# 2. iterate over comparisons
print('formatting info ...')
formatter()

# 3. final print
print('... done.')
