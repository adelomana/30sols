library(pvclust)
library(pheatmap)
library(MASS)
library(gplots)
library(RColorBrewer)
library(entropy)
library(viridis)

# 1. reading data
setwd("~/github/30sols/GREs")
data=read.table('foo.csv',header=FALSE,sep = ",")

# 2. reading and treating metadata
metadata=read.table('classMembership.txt',sep="\t",header=TRUE)
classMembership=data.frame(ClassMembership=metadata$ribosomal.class[1:dim(data)[1]])

working=as.matrix(classMembership)
working[working=='Green, Blue']='black'
working[working=='Green, Magenta']='black'
working[working=='Blue, Magenta']='black'
working[working=='Red, Magenta']='black'

working[working=='Red']='red'
working[working=='Blue']='blue'
working[working=='Green']='green'
working[working=='Magenta']='magenta'
classMembership=as.data.frame(working)

longCoremNames=paste("corem",rownames(data), sep="")
rownames(data)=longCoremNames
rownames(classMembership)=longCoremNames

# computing biclusters
# consider exploring https://bioconductor.org/packages/devel/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html 
#and https://rstudio-pubs-static.s3.amazonaws.com/132435_18fcdc0bd52849f9b61c4b64ced3f22f.html

#permutationIterations=1000

# 2.1. clustering GREs (Vi)
#GREclusterResults=pvclust(data,method.hclust="ward.D2",method.dist="euclidean",nboot=permutationIterations,parallel=TRUE)
#plot(GREclusterResults)
#pvrect(GREclusterResults)
#significantGREclusters=pvpick(GREclusterResults)
#significantGREclusters

# 2.2. clustering corems (i)
#coremClusterResults=pvclust(t(data),method.hclust="ward.D2",method.dist="euclidean",nboot=permutationIterations,parallel=TRUE)
#plot(coremClusterResults)
#pvrect(coremClusterResults)
#significantCoremClusters=pvpick(coremClusterResults)
#significantCoremClusters

# 3. creating figure
#heatmap(as.matrix(data),Rowv=coremClusterResults$hclust$order,Colv=GREclusterResults$hclust$order,col=viridis(256))
#heatmap(as.matrix(data),Rowv=coremClusterResults$hclust$order,Colv=GREclusterResults$hclust$order)

# 4. alternatively, pheatmap
mat_colors <- list(ClassMembership = c('green','black','magenta','blue','red'))
names(mat_colors$ClassMembership) <- unique(classMembership$ClassMembership)
res=pheatmap(data,
             clustering_method='ward.D2',
             annotation_row=classMembership,
             annotation_colors=mat_colors,
             color=viridis(256),
             border_color=NA,
             fontsize          = 8
             #show_colnames=FALSE,
             #show_rownames=FALSE
             )

# consider using entropy as proxy for order
# consider using isa, https://cran.r-project.org/web/packages/isa2/isa2.pdf and https://faculty.washington.edu/dwitten/people.html

