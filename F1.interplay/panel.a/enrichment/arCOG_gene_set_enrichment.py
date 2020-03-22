###
### This script runs gene set enrichment based on arCOG based on annotation from https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/
###

import sys, requests, re
import logging, socket, time

## 1. create a database with DETs, group and arCOG
## 2. run enrichment---Fisher exact test with fdr_bh

def arCOG_mapper():

    
    logger.info('map RS IDs to arCOG')
    arcog2rs = {}

    for arCOG in arCOGs:
        PIs = cog2pi[arCOG]
        RSs = []

        print(arCOG, PIs)
        
        for PI in PIs:
            NP = None; WP = None; RS = None
            try:
                NP = pi2np[PI]
                WP = np2wp[NP]
                RS = wp2rs[WP]
            except:
                pass
            #print(arCOGs.index(arCOG), arCOG, PI, NP, WP, RS)
            if RS is not None:
                RSs.append(RS)
        arcog2rs[arCOG]=RSs

    mapped_RSs = []
    for arCOG in arcog2rs:
        RSs = arcog2rs[arCOG]
        for RS in RSs:
            if RS not in mapped_RSs:
                mapped_RSs.append(RS)
    logger.info('{} RS identifiers have a arCOG mapping'.format(len(mapped_RSs)))

    rs2arcog = {}
    for RS in mapped_RSs:
        rs2arcog[RS] = []
        for arCOG in arcog2rs:
            if RS in arcog2rs[arCOG]:
                rs2arcog[RS].append(arCOG)
        
    return rs2arcog

def arCOG_reader():

    logger.info('map arCOG to protein ID')
    cog2pi = {}
    arCOGs = []
    with open(arCOG_annotations_file,'r') as f:
        for line in f:
            v=line.split(',')
            if len(v) >= 8:
                PI=v[2]
                arCOG=v[6]

                if arCOG not in arCOGs:
                    arCOGs.append(arCOG)

                if arCOG not in cog2pi:
                    cog2pi[arCOG]=[PI]
                else:
                    cog2pi[arCOG].append(PI)
                
    return cog2pi, arCOGs

def wp_fetcher():

    for pi in pi2wp:
        logger.info('fetching info about {}...'.format(pi2wp[pi]))
        
        url = 'https://www.ncbi.nlm.nih.gov/protein/' + pi2wp[pi]
        logger.info(url)
        wpID = None; npID = None
        body = requests.get(url).text.split('\n')
        for line in body:
            if 'identical protein' in line:
                vector = re.split('protein/|\?|Old |\<', line)
                for element in vector:
                    if 'WP_' in element:
                        wpID = element
                    if 'NP_' in element:
                        npID = element
        if (wpID == None) and (npID == None):
                logger.info('problem')
                sys.exit()
        else:
            logger.info('association found\t{}\t{}'.format(wpID, npID))
        
    return None

def np_converter():

    logger.info('map old WP to new WP')
    np2wp = {}
    with open(np_conversion_file, 'r') as f:
        for line in f:
            v = line.split('\t')
            new = v[3]
            old = v[4].replace('\n', '')
            np2wp[old] = new

    return np2wp

def np_reader():

    logger.info('map PI to NP')
    pi2np={}
    with open(PI2NP_file, 'r') as f:
        for line in f:
            v = line.split('|')
            PI = v[1]
            NP = v[3]
            pi2np[PI] = NP
        
    return pi2np

def wp2rs_mapper():

    logger.info('map new WP IDs to RS IDs using GCF')
    wp2rs = {}
    with open(gcf_file, 'r') as f:
        for line in f:
            v = line.split('\t')
            if len(v) > 3:
                info = v[-1]
                if 'VNG_RS' in info and 'WP_' in info:
                    RS = 'VNG_RS' + info.split('VNG_RS')[1].split(';')[0]
                    WP = 'WP_' + info.split('WP_')[1].split(';')[0].split('-')[0]
                    wp2rs[WP] = RS
            
    return wp2rs

#
# MAIN
#

# 0. user defined variables
arCOG_annotations_file='/Volumes/omics4tb2/alomana/projects/TLR/data/arCOG/hsa_nrc1_ar14.arCOG.csv'
PI2NP_file = '/Volumes/omics4tb2/alomana/projects/TLR/data/arCOG/ar14.fa.subset.txt'
np_conversion_file = '/Volumes/omics4tb2/alomana/projects/TLR/data/arCOG/old_new_conversion.txt'
gcf_file = '/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/GCF_000006805.1_ASM680v1_genomic.gff'
DET_dir = '/Users/alomana/github/30sol/F1.interplay/panel.a/results'

# 0.1. logging information
hostname = socket.gethostname()
now = time.strftime('%Y-%m-%d %H:%M:%S')
working_format = '{} {} %(levelname)s | %(message)s'.format(hostname, now)
logging.basicConfig(format=working_format)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

#
# 1. create a database with DETs, group and arCOG

# 1.1. associate protein IDs to arCOG
cog2pi, arCOGs = arCOG_reader()

# 1.2. associate PI to new WP
pi2np = np_reader()

# 1.3. associate old WP to new WP
#wp_fetcher()
np2wp = np_converter()

# 1.4. associate WP to RS
wp2rs = wp2rs_mapper()

# 1.5. associate RS to arCOG
rs2arCOG = arCOG_mapper()

"""
# 1.2. associate geneIDs to color
gene_name2color=DETs_reader()

# 1.3. associate gene names to GI
gene_name2GI=annotation_reader()

# 1.

# run enrichment    
"""
