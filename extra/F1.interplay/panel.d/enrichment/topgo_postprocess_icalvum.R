# Perform post processing of TopGO enrichment to create a table for all groups
post.process <- function(resultsfile=NULL){
  cat("Running Post processing.\n")
  BP <- data.frame()
  CC <- data.frame()
  MF <- data.frame()
  goresults.table <- data.frame()
  
  for(element in names(resultsfile)){
    if(length(resultsfile[[element]]$BP$GO.ID) != 0) {
      BP <- rbind(BP, cbind(group=element,resultsfile[[element]]$BP, GO="BP"))
    }
    
    if(length(resultsfile[[element]]$CC$GO.ID) != 0) {
      CC <- rbind(CC, cbind(group=element,resultsfile[[element]]$CC, GO="CC"))
    }
    
    if(length(resultsfile[[element]]$MF$GO.ID) != 0) {
      MF <- rbind(MF, cbind(group=element, resultsfile[[element]]$MF, GO="MF"))
    }
    
    goresults.table <- rbind(goresults.table, BP, CC, MF)
  }
  goresults.table <- unique(goresults.table)
  return(goresults.table)
}
