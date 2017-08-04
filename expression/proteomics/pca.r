# this is a script that performs PCA on the proteomics samples
library(ggplot2)

setwd("~/30sols/expression/proteomics")

# 1. read the data
lysate.br1.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/lysate.br1.csv',sep=",",header=TRUE,row.names=1)
lysate.br2.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/lysate.br2.csv',sep=",",header=TRUE,row.names=1)
lysate.br3.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/lysate.br3.csv',sep=",",header=TRUE,row.names=1)

enriched.br1.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/rbf.br1.csv',sep=",",header=TRUE,row.names=1)
enriched.br2.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/rbf.br2.csv',sep=",",header=TRUE,row.names=1)
enriched.br3.full=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/rbf.br3.csv',sep=",",header=TRUE,row.names=1)

# 2. treating the data

# 2.1. finding the intersect of names
commonNames=Reduce(intersect,list(rownames(lysate.br1.full),rownames(lysate.br2.full),rownames(lysate.br3.full),rownames(enriched.br1.full),rownames(enriched.br2.full),rownames(enriched.br3.full)))

# 2.2. defining an appropriate data set for PCA
joinedData=data.frame(geneNames=commonNames,
                      lysate.br1.time21=lysate.br1.full[commonNames,"e2.vs..e1_log2_sfc"],
                      lysate.br1.time31=lysate.br1.full[commonNames,"e3.vs..e1_log2_sfc"],
                      lysate.br1.time41=lysate.br1.full[commonNames,"e4.vs..e1_log2_sfc"],
                      #
                      lysate.br2.time21=lysate.br2.full[commonNames,"e2.vs..e1_log2_sfc"],
                      lysate.br2.time31=lysate.br2.full[commonNames,"e3.vs..e1_log2_sfc"],
                      lysate.br2.time41=lysate.br2.full[commonNames,"e4.vs..e1_log2_sfc"],
                      #
                      lysate.br3.time21=lysate.br3.full[commonNames,"e2.vs..e1_log2_sfc"],
                      lysate.br3.time31=lysate.br3.full[commonNames,"e3.vs..e1_log2_sfc"],
                      lysate.br3.time41=lysate.br3.full[commonNames,"e4.vs..e1_log2_sfc"],
                      ###
                      enriched.br1.time21=enriched.br1.full[commonNames,"e2.vs..e1_log2_sfc"],
                      enriched.br1.time31=enriched.br1.full[commonNames,"e3.vs..e1_log2_sfc"],
                      enriched.br1.time41=enriched.br1.full[commonNames,"e4.vs..e1_log2_sfc"],
                      #
                      enriched.br2.time21=enriched.br2.full[commonNames,"e2.vs..e1_log2_sfc"],
                      enriched.br2.time31=enriched.br2.full[commonNames,"e3.vs..e1_log2_sfc"],
                      enriched.br2.time41=enriched.br2.full[commonNames,"e4.vs..e1_log2_sfc"],
                      #
                      enriched.br3.time21=enriched.br3.full[commonNames,"e2.vs..e1_log2_sfc"],
                      enriched.br3.time31=enriched.br3.full[commonNames,"e3.vs..e1_log2_sfc"],
                      enriched.br3.time41=enriched.br3.full[commonNames,"e4.vs..e1_log2_sfc"],
                      #
                      row.names=1)
finalData=t(joinedData)

# 2. perform pca
result=prcomp(finalData,scale=FALSE)
summary(result)
#plot(result$x)
autoplot(result)
