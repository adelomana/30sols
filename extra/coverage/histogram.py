###
### This script builds histograms from the coverage profile text files.
###

import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def figureMaker():

    '''
    This function builds a figure of the coverage of reads over ribo-pt genes.
    '''

    # f.1. read the data

    # f.2. make the figure


    # f.2.3 final figure closing

    #matplotlib.pyplot.xlabel('Time point')
    #matplotlib.pyplot.ylabel('Stoichiometry (log$_2$ ribo-pt)')
    #matplotlib.pyplot.title(condition)

    #matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=3,ncol=2,fontsize=14)

    #figureName='figure.{}.pdf'.format(condition)
    #matplotlib.pyplot.tight_layout()
    #matplotlib.pyplot.savefig(figureName)
    #matplotlib.pyplot.clf()    

    return None

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])

    riboPtNames.sort()
            
    return riboPtNames

###
### MAIN
###

# 0. user defined variables
coverageDir='/Volumes/omics4tb/alomana/projects/TLR/data/coverage/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

timepoints=['tp.1','tp.2','tp.3','tp.4']
replicates=['rep.1','rep.2','rep.3']
experiments=['rbf','trna']

# 1. iterate over ribo-pt genes
riboPtNames=riboPtNamesReader()

for name in riboPtNames:
    
    # 2. build figure
    figureMaker()
