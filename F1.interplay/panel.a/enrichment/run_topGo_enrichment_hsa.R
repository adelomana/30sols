#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("topGO")

library(topGO)
library(multtest)

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
      results <- results[results[,6]<=0.05,] # reomve below significance
      results <- results[results[,4]>1,] # remove below 2 hits
      # correct pvalues
      uncorrected = results[,6]
      corrected = p.adjust(uncorrected, method = 'BH')
      results['adj'] = corrected
      results <- results[results[,7]<=0.05,]
    })
    names(o) <- c("BP", "CC", "MF")
    
    return(o)
  })
}

###
### MAIN
###

# 0. user defined variables
det_folder = '/Users/alomana/github/30sol/F1.interplay/panel.a/results/'
enrichment_folder = '/Users/alomana/github/30sol/F1.interplay/panel.a/enrichment/results/'
labels = c('black.plus', 'black.minus', 'blue', 'dubious', 'green', 'orange', 'red', 'yellow', 'others')
setwd('/Users/alomana/github/30sol/F1.interplay/panel.a/enrichment/')

# 1. iterations
tempo = lapply(labels, function(label){
  input_file = paste(det_folder, 'results.', label, '.txt', sep='')
  output_file = paste(enrichment_folder, 'enrichment.', label, '.GO.DETs.txt', sep='')
  
  print(label)
  print(input_file)
  print(output_file)
  
  # 1.1. read  
  clusters <- read.delim(input_file, sep="\t", header=F)
  mylist <- list()
  mylist[[label]] <- as.vector(clusters$V2)
  
  # 1.2. analysis
  hsa.enrichment <- run.topGO.enrichment(mylist)
  names(hsa.enrichment) <- names(mylist)
  hsa.go.enrichment <- post.process(hsa.enrichment)
  write.table(hsa.go.enrichment, file=output_file, sep="\t", row.names=F)
  
  print('----')
})