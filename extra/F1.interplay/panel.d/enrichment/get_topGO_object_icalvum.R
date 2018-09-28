# Get  topGO data Object  [6]    
get_topGO_object <- function(genes,gene2go,ontology=c("BP","MF","CC")) {
  require(topGO)
  cat("Getting topGO object","\n")
  # genes is a vector containing genes of interest
  geneList <- factor(as.integer(names(gene2go)%in%genes))
  names(geneList) <- names(gene2go)
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go)
  #GOdata can be used directly for analysis, e.g. 
  test <- runTest(GOdata,algorithm="classic",statistic="fisher")
  results <- GenTable(GOdata,test,topNodes=20)
  return(GOdata)
}