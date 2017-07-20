library(pvclust)
library(pheatmap)
library(MASS)
library(gplots)
library(RColorBrewer)
library(entropy)

# 1. reading data
setwd("~/github/30sols/GREs")
data=read.table('foo.csv',header=FALSE,sep = ",")

metadata=read.table('classMembership.txt',sep="\t",header=TRUE)
classMembership=data.frame(ClassMembership = metadata$ribosomal.class[1:dim(data)[1]])
longCoremNames=paste("corem",rownames(data), sep="")
rownames(data)=longCoremNames
rownames(classMembership)=longCoremNames

# computing biclusters
# consider exploring https://bioconductor.org/packages/devel/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html and https://rstudio-pubs-static.s3.amazonaws.com/132435_18fcdc0bd52849f9b61c4b64ced3f22f.html

permutationIterations=1000

# 2.1. clustering GREs (Vi)
GREclusterResults=pvclust(data,method.hclust="ward.D2",method.dist="euclidean",nboot=permutationIterations,parallel=TRUE)
plot(GREclusterResults)
pvrect(GREclusterResults)
significantGREclusters=pvpick(GREclusterResults)
significantGREclusters

# 2.2. clustering corems (i)
coremClusterResults=pvclust(t(data),method.hclust="ward.D2",method.dist="euclidean",nboot=permutationIterations,parallel=TRUE)
plot(coremClusterResults)
pvrect(coremClusterResults)
significantCoremClusters=pvpick(coremClusterResults)
significantCoremClusters

# 3. creating figure
heatmap(as.matrix(data),Rowv=coremClusterResults$hclust$order,Colv=GREclusterResults$hclust$order,col=viridis(256))

# 4. alternatively, pheatmap
res=pheatmap(data,clustering_method='complete',annotation_row=classMembership)

# consider using entropy as proxy for order
# consider using isa, https://cran.r-project.org/web/packages/isa2/isa2.pdf and https://faculty.washington.edu/dwitten/people.html

