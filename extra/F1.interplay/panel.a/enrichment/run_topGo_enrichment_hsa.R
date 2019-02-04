clusters <- read.delim("/Users/alomana/github/30sol/extra/F1.interplay/panel.d/results/results.black.minus.txt", sep="\t", header=F)
mylist <- list()
mylist[["blackMinus"]] <- as.vector(clusters$V1)

# enrichment analysis
run.topGO.enrichment <-function(my.members){
  source("~/scratch/get_topGO_object_hsa.R")
  source("~/scratch/gene2go_function_hsa.R")
  source("~/scratch/topgo_postprocess_hsa.R")
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
      test <- runTest(tmp,algorithm="classic",statistic="fisher")
      results <- GenTable(tmp,test,topNodes=length(test@score))
      results <- results[results[,6]<=0.05,]
    })
    names(o) <- c("BP", "CC", "MF")
    
    return(o)
  })
}

library(topGO)

hsa.enrichment <- run.topGO.enrichment(mylist)
names(hsa.enrichment) <- names(mylist)
hsa.go.enrichment <- post.process(hsa.enrichment)
write.table(hsa.go.enrichment, file="~/scratch/black.minus.DET_GO_Enrichment.txt", sep="\t", row.names=F)
