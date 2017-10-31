import os,sys

def genomeIndexer():

    '''
    this function creates the genome index. Should be run only once.
    '''

    flag1=' --runMode genomeGenerate'
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --genomeDir %s'%genomeIndexDir
    flag4=' --genomeFastaFiles %s'%genomeFastaFile
    flag5=' --sjdbGTFfile %s'%genomeAnnotationFile
    flag6=' --sjdbGTFfeatureExon gene --sjdbGTFtagExonParentTranscript locus_tag --sjdbGTFtagExonParentGene Name --sjdbOverhang 30'
    flag7=' --genomeSAindexNbases 10'

    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag7+flag5+flag6

    print
    print cmd
    print
    os.system(cmd)

    return None

# 0. defining several input/output paths
STARexecutable='/Users/arjun/STAR/bin/MacOSX_x86_64/STAR'
genomeIndexDir='/Users/arjun/data/genomeIndex_noAnn'
genomeFastaFile='/Users/arjun/data/ncbi/GCF_000006805.1_ASM680v1_genomic.fna'
genomeAnnotationFile='/Users/arjun/data/ncbi/halo.genomic.archive.gff'
#smallgenome=10
numberOfThreads=4

# 2. making genome indexes
print 'making genome index...'
genomeIndexer()
