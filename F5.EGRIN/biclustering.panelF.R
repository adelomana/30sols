#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pheatmap")
#BiocManager::install("viridis")

library(pheatmap)
library(viridis)

here=getwd()
setwd(here)

# 1. reading and treating data
data=as.matrix(read.table('/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.clean.csv',header=TRUE,sep = ",",row.names=1))
valid_indices = c()
for  (index in c(1:dim(data)[2]))
{
  sumGREs = sum(data[, index])
  if (sumGREs > 100){
    valid_indices = c(valid_indices, index)
  }
}

data = data[, valid_indices]
head(data)
data = log10(data+1)
head(data)

# 2. reading and treating metadata
metadata=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/annotation/coremClasses.csv',sep=",",header=TRUE,row.names=1)
plottingCorems=rownames(data)
selectedMetadata=data.frame(corems=plottingCorems,ribo.class=metadata[plottingCorems,],row.names=1)

annotationColors=list(ribo.class=levels(selectedMetadata$ribo.class))
names(annotationColors$ribo.class)=levels(selectedMetadata$ribo.class)

# 3. making a heatmap
res=pheatmap(data,clustering_method='ward.D2',show_rownames = TRUE,show_colnames = TRUE,color=viridis(256),fontsize=8,border_color=NA,annotation_row=selectedMetadata,annotation_colors=annotationColors)

# # 4. GREs
# GREs are directly grabbed from the web, e.g., http://egrin2.systemsbiology.net/gres/64091/hal_32