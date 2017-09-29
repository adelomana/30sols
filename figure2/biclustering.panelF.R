library(pheatmap)
library(viridis)

# 1. reading and treating data
setwd("~/github/30sols/GREs")
data=as.matrix(read.table('/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.clean.csv',header=TRUE,sep = ",",row.names=1))
data[data>50]=50
dim(data)

# 2. reading and treating metadata
metadata=read.table('/Volumes/omics4tb/alomana/projects/TLR/data/annotation/coremClasses.csv',sep=",",header=TRUE,row.names=1)
plottingCorems=rownames(data)
selectedMetadata=data.frame(corems=plottingCorems,ribo.class=metadata[plottingCorems,],row.names=1)

annotationColors=list(ribo.class=levels(selectedMetadata$ribo.class))
names(annotationColors$ribo.class)=levels(selectedMetadata$ribo.class)

# 3. making a heatmap
res=pheatmap(data,clustering_method='ward.D2',show_rownames = TRUE,show_colnames = TRUE,color=viridis(256),fontsize=8,border_color=NA,annotation_row=selectedMetadata,annotation_colors=annotationColors)

# 4. making pie charts
slices=c(1)
lbls=c('red')
pie(slices,labels=lbls,main='pie 1',col=lbls)

slices=c(7,9)
lbls=c('magenta','blue')
pie(slices,labels=lbls,main='pie 2',col=lbls)

slices=c(9,4,10,1)
lbls=c('red','green','magenta','blue')
pie(slices,labels=lbls,main='pie 3',col=lbls)

slices=c(9,1)
lbls=c('magenta','blue')
pie(slices,labels=lbls,main='pie 4',col=lbls)