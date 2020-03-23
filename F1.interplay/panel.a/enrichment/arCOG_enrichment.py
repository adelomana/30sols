###
### This script performs functional enrichment based on arCOG
###

import sys
import logging, socket, time
import scipy, scipy.stats
import statsmodels, statsmodels.stats, statsmodels.stats.multitest

def annotation_reader():

    annotation = {}
    with open(functional_annotation_file, 'r') as f:
        for line in f:
            v = line.split()
            rs = v[0]
            cog = v[1]

            if rs not in annotation:
                annotation[rs] = [cog]
            else:
                annotation[rs].append(cog)

    # summary
    retrieved_genes_rank = len(list(annotation.keys()))
    cogs = []
    for rs in annotation:
        values = annotation[rs]
        for value in values:
            if value not in cogs:
                cogs.append(value)
    retrieved_cogs_rank = len(cogs)
    logger.info('retrieved RS: {}; retrieved arCOGs: {}'.format(retrieved_genes_rank, retrieved_cogs_rank))

    return annotation, retrieved_genes_rank, retrieved_cogs_rank

def DET_reader(color):

    DETs = []; subset = []
    working_file = '{}results.{}.txt'.format(DET_folder, color)
    print(working_file)
    with open(working_file, 'r') as f:
        for line in f:
            v = line.split()
            rs = v[2]
            rs = rs.replace('gene-', '')
            rs = rs.replace('VNG', 'VNG_')
            if rs not in DETs:
                DETs.append(rs)
                if rs in annotation:
                    subset.append(rs)
    logger.info('retrieved DETs: {}; with annotation: {}'.format(len(DETs), len(subset)))
    
    return subset

def enrichment(color):

    # f.1. define categories in group
    logger.info('defining functional categories in color {}'.format(color))
    categories = []
    for DET in DETs:
        if DET in annotation:
            for element in annotation[DET]:
                if element not in categories:
                    categories.append(element)
    logger.info('detected {} categories in color {}'.format(len(categories), color))

    # f.2. Fisher exact test
    #        setA  setB
    # hit
    # no hit
    #
    # hit in setA, hit in setB, no hit in setA, no hit in setB

    logger.info('computing Fisher exact test for each category')

    pvals = []
    
    for category in categories:

        
        # define elements of that category in the color
        hits_in_a = 0; opp1 = 0
        print(len(DETs))
        for DET in DETs:
            if category in annotation[DET]:
                hits_in_a = hits_in_a + 1
            else:
                opp1 = opp1 + 1
        print(hits_in_a, opp1)
                

        # define elements of that category in the rest of the genome that is not the color
        hits_in_B = 2

        # define elements in color that do not have that category
        no_hits_in_A = 2

        # define elements in the rest of the genome that is not the color, that do not have that category
        no_hits_in_B = 5

        # perform test
        hits_in_a = 8
        oddsratio, pvalue = scipy.stats.fisher_exact([[hits_in_a, hits_in_b], [no_hits_in_A, no_hits_in_b]])
        print(oddsratio, pvalue)
        pvals.append(pvalue)

    # f.3. multiple testing correction
    logger.info('correcting for multiple-testing')
    reject, pvals_corrected = statsmodels.stats.multitest.multipletests(pvals, alpha=0.1, method='fdr_bh')
    print(pvals)
    print(pvals_corrected)
    print(reject)
    sys.exit()
    
    return None

###
### MAIN
###

# 0. user-defined variables
functional_annotation_file = '/Volumes/omics4tb2/alomana/projects/TLR/results/arCOG/arCOG.final.annotation.txt'
DET_folder = '/Users/alomana/github/30sol/F1.interplay/panel.a/results/'
colors = ['blue', 'orange']

# 0.1. logging information
hostname = socket.gethostname()
now = time.strftime('%Y-%m-%d %H:%M:%S')
working_format = '{} {} %(levelname)s | %(message)s'.format(hostname, now)
logging.basicConfig(format=working_format)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# 1. read information
logger.info('reading files')
annotation, retrieved_genes_rank, retrieved_cogs_rank = annotation_reader()
print(annotation)
# 2. analysis
logger.info('analysis')

# 2.1. iterate over colors
for color in colors:
    DETs = DET_reader(color)
    enrichment(color)
