###
### This script computes normalized counts for all samples using DESeq2.
###

library('DESeq2')
library('stringr')

# 0. user defined variables
setwd("~/github/aukera/extra/expression/pipeline.star.htseq-count.deseq2")
countsDir="/Volumes/omics4tb/alomana/projects/TLR/data/counts"
DESeqNormalizedCountsFile=paste('/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv')

# 1. handle data reading
sampleFiles=dir(file.path(countsDir))

# 1.1. selecting data files and data names
sampleNames=str_sub(sampleFiles,1,-5)

# 1.2. defining a random design
timeCondition=sapply(strsplit(sampleNames,split='.',fixed=TRUE),function(x) (paste('tp',x[5],sep='.')))

# 2. formatting variables
sampleTable=data.frame(sampleName=sampleNames,fileName=sampleFiles,condition=timeCondition)
dds=DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=countsDir,design= ~ condition)

# 3. filtering low abundance transcripts
keep=rowSums(counts(dds)) >= 10
dds=dds[keep,]

# 4. normalized values
vsd=vst(dds,blind=FALSE)
rld=rlog(dds,blind=FALSE)

# 5. exporting normalized counts
write.csv(as.data.frame(assay(rld)),file=DESeqNormalizedCountsFile,quote=FALSE)
