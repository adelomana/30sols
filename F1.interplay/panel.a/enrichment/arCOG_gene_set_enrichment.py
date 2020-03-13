###
### This script runs gene set enrichment based on arCOG
###

## 1. create a database with DETs, group and arCOG
## 2. run enrichment---Fisher exact test with fdr_bh

def arCOG_reader():

    '''
    Read arCOG annotations.
    '''

    GI2COG={}
    with open(arCOG_annotations_file,'r') as f:
        for line in f:
            v=line.split(',')
            if len(v) >= 8:
                GI=v[2]
                arCOG=v[6]

                if GI not in GI2COG:
                    GI2COG[GI]=[arCOG]
                else:
                    GI2COG[GI].append(arCOG)
                
    return GI2COG

#
# MAIN
#

# 0. user defined variables
arCOG_annotations_file='/Volumes/omics4tb2/alomana/projects/TLR/data/arCOG/hsa_nrc1_ar14.arCOG.csv'
DETs_file=''
annotation_file='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.1.cs.NC_002608.1.fasta'

# 1. create a database with DETs, group and arCOG

# 1.1. associate protein IDs to arCOG
GI2COG=arCOG_reader()

# 1.2. associate geneIDs to color
gene_name2color=DETs_reader()

# 1.3. associate gene names to GI
gene_name2GI=annotation_reader()

# run enrichment    
