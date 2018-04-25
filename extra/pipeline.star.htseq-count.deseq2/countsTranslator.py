###
### This script transforms the new annotation to the old annotation in the counts file.
###


def synonymReader():

    '''
    This function builds a dictionary for conversion between old locus names and current locus names.
    '''

    synonyms={}
    with open(genomeAnnotationFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1]
                currentLocus=info.split('Parent=')[1].split(';')[0]
                print(currentLocus)

    return synonyms

###
### MAIN
###

# 0. user defined variables
genomeAnnotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

# 1. reading data
synonyms=synonymReader()

# 2. performing the analysis
