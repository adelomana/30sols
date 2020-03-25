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

    annotation_names = {}
    with open(arCOG_names, 'r') as f:
        for line in f:
            vector = line.split('\t')
            ID = vector[0]
            info = vector[3]
            annotation_names[ID] = info
    annotation_names['arCOG_null'] = 'Unknown function'

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

    return annotation, annotation_names, retrieved_genes_rank, retrieved_cogs_rank

def DET_reader(color):

    DETs = []; subset = []
    working_file = '{}results.{}.txt'.format(DET_folder, color)
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
                #else:
                #    print('{} not found'.format(rs))
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

    pvals = []; hits = []; expected = []
    
    for category in categories:

        # define elements of a functional category in the color
        hits_in_a = 0; opp1 = 0
        #print('lenDETs', len(DETs))
        for DET in DETs:
            if category in annotation[DET]:
                hits_in_a = hits_in_a + 1
            else:
                opp1 = opp1 + 1
        #print('hits_ina', hits_in_a, opp1)

        # define elements of that functional category in the rest of the genome that is not the color
        hits_in_b = 0
        for rs in annotation:
            if category in annotation[rs]:
                if rs not in DETs:
                    hits_in_b = hits_in_b + 1
        #print('hits in bkground', hits_in_b)

        # define elements in color that do not have that category
        no_hits_in_a = len(DETs) - hits_in_a
        #print('nohitsinA', no_hits_in_a)

        # define elements in the rest of the genome that is not the color, that do not have that category
        no_hits_in_b = 0
        for rs in annotation:
            if rs not in DETs:
                if category not in annotation[rs]:
                    no_hits_in_b = no_hits_in_b + 1
        #print('nohitsinB', no_hits_in_b)
        #print('total genome', len(list(annotation.keys())))

        # perform test
        oddsratio, pvalue = scipy.stats.fisher_exact([[hits_in_a, hits_in_b], [no_hits_in_a, no_hits_in_b]])
        #print(category, hits_in_a, hits_in_b, no_hits_in_a, no_hits_in_b, pvalue)
        
        pvals.append(pvalue)
        hits.append(hits_in_a)
        expect = ((hits_in_a + hits_in_b) / len(list(annotation.keys()))) * len(DETs)
        expected.append(expect)

    # f.3. multiple testing correction
    logger.info('correcting for multiple-testing')
    output = statsmodels.stats.multitest.multipletests(pvals, alpha=0.1, method='fdr_bh')
    alt_hypotheses = output[0]
    pvals_corrected = output[1]

    # f.4. store results
    set_size = len(DETs)
    enrichment_file = results_dir + color + '.txt'
    with open(enrichment_file, 'w') as f:
        for i in range(len(categories)):
            if (alt_hypotheses[i] == True) and (pvals[i] < 0.05) or (categories[i] == 'arCOG_null'):
                info = annotation_names[categories[i]]
                hit_rank = hits[i]
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(color, set_size, categories[i], info, hit_rank, pvals[i], pvals_corrected[i]))
    
    return None

###
### MAIN
###

# 0. user-defined variables
functional_annotation_file = '/Volumes/omics4tb2/alomana/projects/TLR/results/arCOG/arCOG.final.annotation.txt'
DET_folder = '/Users/alomana/github/30sol/F1.interplay/panel.a/results/'
colors = ['black.plus', 'black.minus', 'blue', 'dubious', 'green', 'orange', 'red', 'yellow']
results_dir = '/Volumes/omics4tb2/alomana/projects/TLR/results/arCOG/enrichment/'
arCOG_names = '/Volumes/omics4tb2/alomana/projects/TLR/data/arCOG/ar14.arCOGdef19.tab'

# 0.1. logging information
hostname = socket.gethostname()
now = time.strftime('%Y-%m-%d %H:%M:%S')
working_format = '{} {} %(levelname)s | %(message)s'.format(hostname, now)
logging.basicConfig(format=working_format)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# 1. read information
logger.info('reading files')
annotation, annotation_names, retrieved_genes_rank, retrieved_cogs_rank = annotation_reader()

# 2. analysis
logger.info('analysis')

# 2.1. iterate over colors
for color in colors:
    DETs = DET_reader(color)
    enrichment(color)
