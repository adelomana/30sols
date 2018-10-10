###
### This script defines differentially expressed transcripts using DESeq2 and saves normalized counts for only those samples.
###

library('DESeq2')
library('stringr')

# 0. user defined variables
tag='rbf' # trna / rbf
earlyTimepointName='tp.1'
lateTimepointName='tp.2'
comparison='condition_tp.2_vs_tp.1'

setwd("~/scratch/")
countsDir="/Volumes/omics4tb/alomana/projects/TLR/data/counts"

DESeqResultsFile=paste("/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/significance",tag,comparison,"csv",sep='.')
DESeqNormalizedCountsFile=paste('/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts',tag,comparison,'csv',sep='.')

# 1. handle data reading
sampleIDs=dir(file.path(countsDir))

# 1.1. selecting data files and data names
sampleSubset=sampleIDs[str_detect(sampleIDs,tag)]
ti=sampleSubset[str_detect(sampleSubset,earlyTimepointName)]
tf=sampleSubset[str_detect(sampleSubset,lateTimepointName)]
samples=c(ti,tf)
timeCondition=sapply(strsplit(samples,split='.',fixed=TRUE),function(x) (paste('tp',x[5],sep='.')))

# 2. formatting variables
sampleTable=data.frame(sampleName=samples,fileName=samples,condition=timeCondition)
dds=DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=countsDir,design= ~ condition)

# 3. filtering low abundance transcripts
keep=rowSums(counts(dds)) >= 10
dds=dds[keep,]

# 4. performing differential expression
dds=DESeq(dds)
res=results(dds)
summary(res)

resultsNames(dds)
resLFC=lfcShrink(dds,coef=comparison)

res05=results(dds,alpha=0.05)
summary(res05)

# 5. plotting
plotMA(res)
plotMA(resLFC)

plotCounts(dds, gene=which.min(res$padj),intgroup="condition")
plotCounts(dds, gene=which.max(res$padj),intgroup="condition")

# 6. saving CSV files
resOrdered=res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),file=DESeqResultsFile)

# 7. exporting normalized counts
#values=counts(dds,normalized=TRUE)
#write.csv(as.data.frame(values),file=DESeqNormalizedCountsFile)
