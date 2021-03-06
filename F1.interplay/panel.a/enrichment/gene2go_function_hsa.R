gene2go.object <-function(file=NULL){
  file = "/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/genome.annotation.txt"
  cat("Loading genome annotation file from ", file, "\n")
  Proteome <- read.delim(file, header=T, sep="\t")
  cat("creating gene2go file ", "\n")
  
  gene2go = list()
  print(summary(length(Proteome$sysName)))
  for(gene in Proteome$sysName){
    name = sub("_", "", gene)
    if(length(strsplit(as.character(Proteome[which(Proteome$sysName== gene),"GO"]), split = ",", fixed=T)[[1]]) > 0){
      gene2go[gene] = strsplit(as.character(Proteome[which(Proteome$sysName== gene),"GO"]), split = ",", fixed=T)
    }
  }
  return(gene2go)
}