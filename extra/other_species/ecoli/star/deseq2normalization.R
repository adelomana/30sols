###
### This script computes normalized counts for all samples using DESeq2. Normalization is carried out for RNA-seq and Ribo-seq separately.
###

# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library('DESeq2')
library('stringr')

# 0. user defined variables
tag='rbf' # trna / rbf

setwd("~/scratch/")
countsDir="/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/counts"
DESeqNormalizedCountsFile=paste('//Volumes/omics4tb2/alomana/projects/TLR/data/ecoli//DESeq2/normalizedCounts',tag,'csv',sep='.')

# 1. handle data reading
sampleFiles=dir(file.path(countsDir))
sampleSubset=sampleFiles[str_detect(sampleFiles,tag)]

# 1.1. selecting data files and data names
sampleNames=str_sub(sampleSubset,1,-5)

# 1.2. defining a random design
timeCondition=sapply(strsplit(sampleNames,split='.',fixed=TRUE),function(x) (paste('tp',x[5],sep='.')))

# 2. formatting variables
sampleTable=data.frame(sampleName=sampleNames,fileName=sampleSubset,condition=timeCondition)
dds=DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=countsDir,design= ~ condition)

# 3. filtering low abundance transcripts
keep=rowSums(counts(dds)) >= 10
dds=dds[keep,]

# 4. normalized values
vsd=vst(dds,blind=FALSE)
rld=rlog(dds,blind=FALSE)

# 5. exporting normalized counts
write.csv(as.data.frame(assay(rld)),file=DESeqNormalizedCountsFile,quote=FALSE)
