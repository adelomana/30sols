### script to generate Figures 1A and B of manuscript

bootplot <- function(genelist, network = 'corem', transpose = TRUE) {
  library(pvclust)
  library(RColorBrewer)
  corem.list <- lapply(genelist, paste(network,'.finder', sep=''))
  unique.corems <- unique(unlist(corem.list))
  heat <- matrix(nrow=length(genelist), ncol=length(unique.corems), dimnames=list(riboSubunits[genelist,3], unique.corems))
  for(i in 1:length(corem.list)) {heat[i,] <- unique.corems %in% corem.list[[i]]}
  heatR <- heat + 0
  heatR <- heatR[rowSums(heatR) != 0,]
  if(transpose == TRUE) {
    heat <- t(heatR)
  } 
  result <- pvclust(heat, method.dist = 'binary', method.hclust = 'average', nboot = 10000)
  plot(result, print.pv = FALSE, print.num = FALSE, font = 1)
  pvrect(result, alpha = 0.95, border = 8, lwd = 2)
  heatmap(heatR, xlab="Corem", ylab="protein", col=brewer.pal(3, 'Greens'), scale="none", distfun = function(x) dist(x, method = 'binary'), hclustfun = function(x) hclust(x, method = 'average'))
}
