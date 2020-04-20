# alorenzetti 2020

# this script will take the network file
# and perform the functional categorization using the provided annotation file

# setting wd watch out
setwd("/Users/alan/gdrive/documentos/doutorado/isb/20200113-cytoscape/")

# loading libs
library(pacman)
libs = c("rtracklayer",
          "tidyverse",
          "openxlsx")
p_load(char=libs)

# # loading list of ribosomal protein names (provided by Adrián)
# riboNames= read_delim("ribosomalGeneNames.txt", delim = "\t", col_names = T)
# colnames(riboNames) = gsub(" ", "_", colnames(riboNames))

# dataset3 from https://www.pnas.org/highwire/filestream/592990/field_highwire_adjunct_files/7/11663Dataset_3.sif.txt
nodes = read_delim("Dataset_3.sif", delim = "\t", col_names = F)

# getting the list of falsePositives pulled down by ProteinA
# (when proteinA was used as a solo bait)
falsepos = nodes[which(nodes[,1] == "ProteinA"), "X3"]
falsepos = falsepos$X3 %>% as.character()
falsepos = c(falsepos, "ProteinA")

# if the bait is proteinA, remove from dataset
# but before doing it, let s store these proteins
pAlist = nodes[which(nodes[,1] == "ProteinA"),]
colnames(pAlist)[3] = "old_locus_tag"

# now removing proteins baited by proteinA
nodes = nodes[-which(nodes[,1] == "ProteinA"),]

# if the bait pullsdown a falsepositive or proteinA, remove the entry from dataset
# but before doing that, lets make a backup (will be used afterwards)
nodesbckp = nodes
nodes = nodes[!nodes$X3 %in% falsepos,]

# write new parsed file
write_delim(nodes, "Dataset_3_noFalsePositives.sif", delim="\t", col_names = F)

# getting the list of all nodes so we can classify them
allnodes = c(nodes$X1, nodes$X3) %>% unique()
allnodes = enframe(allnodes) %>% select(value)
colnames(allnodes) = "old_locus_tag"
df = allnodes

# loading NCBI old annotation file
ncbi = rtracklayer::import("~/gdrive/dlsm/misc/Hsalinarum.gff")

ncbiCDS = subset(ncbi, type == "CDS")
ncbiCDS = as_tibble(ncbiCDS)
ncbiCDS = ncbiCDS %>% dplyr::select(Parent, product)
ncbiCDS$Parent = ncbiCDS$Parent %>% as.character()

ncbigenes = subset(ncbi, type == "gene")
ncbigenes = as_tibble(ncbigenes)
ncbigenes = ncbigenes %>% dplyr::select(ID, locus_tag, old_locus_tag)
ncbigenes$old_locus_tag = ncbigenes$old_locus_tag %>% sub("^.*,", "", .) %>% sub("VNG_","VNG",.)

ncbiCDS = left_join(ncbiCDS, ncbigenes, by = c("Parent" = "ID"))

# loading third party gene annotation file
tpa = rtracklayer::import("~/gdrive/dlsm/de_analysis/misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
CDS = subset(tpa, type == "CDS")
CDS = as_tibble(CDS)
CDS = CDS %>% dplyr::select(locus_tag, product)

genes = subset(tpa, type == "gene")
genes = as_tibble(genes)
genes = genes %>% dplyr::select(locus_tag)

CDS = left_join(CDS, genes, by = "locus_tag")
CDS$locus_tag = CDS$locus_tag %>% sub("_", "", .)

# full joining both dataframes (all non-matching elements will be kept)
join = full_join(ncbiCDS, CDS, by = c("old_locus_tag" = "locus_tag"))

# loading kegg and uniprot additional info
# getting kegg info
if(file.exists("keggSetHalDf.RData")){
  load("keggSetHalDf.RData")

} else {
  keggSetHal = kegg.gsets(species = "hal", id.type = "kegg")
  keggSetHal = keggSetHal$kg.sets

  keggSetHalDf = keggSetHal %>%
    unlist %>%
    tibble::enframe()

  keggSetHalDf$name = sub(keggSetHalDf$name, pattern = "[0-9]{1,}$", replacement = "")
  colnames(keggSetHalDf) = c("KEGGpathway", "locus_tag")

  keggSetHalDf %<>%
    dplyr::group_by(locus_tag) %>%
    dplyr::summarise(KEGGpathway = base::paste(KEGGpathway, collapse = "; "))

  save(keggSetHalDf, file="keggSetHalDf.RData")
}

# getting uniprot info
# uniprot is not letting us get all the queried terms
# we have to split the queries by 100 or 150 each time
if(file.exists("uniprotHal.RData")){
  load("uniprotHal.RData")

} else {
  uniprotHal = UniProt.ws(taxId=64091)
  columns = c("UNIPROTKB", "GO")

  # this only has to be done when using the custom
  # annotation by pfeiffer et al. 2019
  CDSfilt = CDS[!paste0("hal:", CDS$locus_tag) %>% grepl("a$", .),]
  CDSfilt = CDS[!paste0("hal:", CDS$locus_tag) %>% grepl("b$", .),]

  m = CDSfilt %>% dim %>% .[1]
  v = seq(from=1, to=m, by=100)
  if(v[v %>% length]){v[v %>% length + 1] = m}
  table=NULL

  for(i in seq(1, length(v)-1)){
    if(i == 1){
      chunk = UniProt.ws::select(uniprotHal,
                                 keys=paste0("hal:", CDSfilt$locus_tag)[v[i]:v[i+1]],
                                 columns = columns,
                                 keytype = "KEGG")

      table = bind_rows(table,chunk)

    }else{
      chunk = UniProt.ws::select(uniprotHal,
                                 keys=paste0("hal:", CDSfilt$locus_tag)[(v[i]+1):(v[i+1])],
                                 columns = columns,
                                 keytype = "KEGG")

      table = bind_rows(table,chunk)
    }
  }

  uniprotHal = table
  uniprotHal = uniprotHal[!is.na(uniprotHal$KEGG),]
  uniprotHal = uniprotHal[uniprotHal$KEGG %>% duplicated == FALSE,]
  uniprotHal$KEGG = sub("^hal:", "", uniprotHal$KEGG)

  save(uniprotHal, file="uniprotHal.RData")
}

# removing underscore from locus_tag of kegg and uniprot datasets
keggSetHalDf$locus_tag = keggSetHalDf$locus_tag %>% sub("VNG_", "VNG", .)
uniprotHal$KEGG = uniprotHal$KEGG %>% sub("VNG_", "VNG", .)

# funcat based on old_locus_tag
allnodes = left_join(allnodes, join, by = c("old_locus_tag" = "old_locus_tag")) %>%
  select(old_locus_tag, product.x, product.y)
colnames(allnodes)[2:3] = c("productNCBIRefSeq", "productPfeiffer2019")

allnodes = left_join(allnodes, keggSetHalDf, by = c("old_locus_tag" = "locus_tag"))
allnodes = left_join(allnodes, uniprotHal, by = c("old_locus_tag" = "KEGG"))

# saving results
write_delim(allnodes, "resultsFuncat.tsv", delim="\t", col_names = T)

# doing the classification for the proteins baited by proteinA
# (false pos)
pAlistClas = left_join(pAlist, join, by = c("old_locus_tag" = "old_locus_tag")) %>%
  select(old_locus_tag, product.x, product.y)
colnames(pAlistClas)[2:3] = c("productNCBIRefSeq", "productPfeiffer2019")

# using the pAlistClas object above, I was able to see that
# VNG2657G (S7) and VNG1715G (S5) were baited by proteinA
# reasoning that those are baited by proteinA, and they are
# bound to "third" ribosome proteins
# let's check which proteins baited S7 and S5

# reading tf labels
tfLabels = read_csv("tfLabels.csv", col_names = T)

bind_rows(nodesbckp[nodesbckp$X3 == "VNG2657G",] %>%
  left_join(.,tfLabels, by = c("X1" = "name")),
nodesbckp[nodesbckp$X3 == "VNG1715G",] %>%
  left_join(.,tfLabels, by = c("X1" = "name"))) %>%
  pull("label") %>% sort() %>% unique()
