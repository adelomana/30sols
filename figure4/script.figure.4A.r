                                        # this is a script that performs PCA on the transcriptomics samples (mRNA and footprints)

data <- plotPCA(rld_new, intgroup = c( "sampleType", "timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")

ggplot(data, aes(PC1, PC2, color=sampleType, shape=timepoint)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
