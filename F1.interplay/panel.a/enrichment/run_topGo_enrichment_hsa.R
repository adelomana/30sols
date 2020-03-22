#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("topGO")

library(topGO)

# enrichment analysis function
run.topGO.enrichment <-function(my.members){
  source("get_topGO_object_hsa.R")
  source("gene2go_function_hsa.R")
  source("topgo_postprocess_hsa.R")
  gene2go <- gene2go.object()
  
  # Enrichment analysis 
  bc1.go <- lapply(seq(1,length(my.members)), function(first){
    # do not test groups with less than 2 members
    if(length(my.members[[first]])<2) {    ## added DJR
      cat("\n Not analyzing group ", first, "... \n")   ## added DJR
      return(NULL) ## added DJR
    }   ## added DJR
    # skip bicluster with less than two genes mapping to go terms
    if(sum(is.element(my.members[[first]], names(gene2go))) <2){
      cat("\n Not analyzing group ", first, "... \n")
      return(NULL)
    }
    
    # Do enrichment for all GO categories for all biclusters
    o <- lapply(c("BP", "CC", "MF"), function(second){
      cat("\n Now analyzing group ", first, "with category ", second, "... \n")
      tmp = get_topGO_object(my.members[[first]], gene2go, second)
      test <- runTest(tmp,algorithm="classic",statistic="fisher") # no multiple testing yet!!!
      results <- GenTable(tmp,test,topNodes=length(test@score))
      results <- results[results[,6]<=0.05,]
    })
    names(o) <- c("BP", "CC", "MF")
    
    return(o)
  })
}

###
### MAIN
###

# 0. user defined variables
input_file = '/Users/alomana/github/30sol/F1.interplay/panel.a/results/results.orange.txt'
output_file = '/Users/alomana/github/30sol/F1.interplay/panel.a/enrichment/results/enrichment.orange.GO.DETs.txt'
label = 'black_minus'

# 1. reading
clusters <- read.delim(input_file, sep="\t", header=F)
mylist <- list()
mylist[[label]] <- as.vector(clusters$V2)

# 2. analysis
hsa.enrichment <- run.topGO.enrichment(mylist)
names(hsa.enrichment) <- names(mylist)
hsa.go.enrichment <- post.process(hsa.enrichment)
write.table(hsa.go.enrichment, file=output_file, sep="\t", row.names=F)
