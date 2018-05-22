###
### this script defines differentially expressed transcripts using Sleuth
###

library('sleuth')
library('stringr')

# 0. user defined variables
setwd("~/github/aukera/extra/expression/pipeline.kallisto")
resultsDir='/Volumes/omics4tb/alomana/projects/TLR/data/kallisto1e3'
sleuthResultsRBF='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/sleuthResultsRBF.csv'
sleuthResultsRNA='/Volumes/omics4tb/alomana/projects/TLR/data/sleuth1e3/sleuthResultsRNA.csv'
metadataFileRBF='sleuth.metadata.rbf.csv'
metadataFileRNA='sleuth.metadata.rna.csv'

# 1. handle data reading
sampleIDs=dir(file.path(resultsDir))
kallistoDirs=file.path(resultsDir,sampleIDs)

# 1.1. working with RBF
sampleSubset=sampleIDs[str_detect(sampleIDs,"rbf")]
t1=sampleSubset[str_detect(sampleSubset,"tp.1")]
t4=sampleSubset[str_detect(sampleSubset,"tp.4")]
sampleIDsRBF=c(t1,t4)

dirsSubset=kallistoDirs[str_detect(kallistoDirs,"rbf")]
t1=dirsSubset[str_detect(dirsSubset,"tp.1")]
t4=dirsSubset[str_detect(dirsSubset,"tp.4")]
kallistoDirsRBF=c(t1,t4)

# 1.2. working with RNA
sampleSubset=sampleIDs[str_detect(sampleIDs,"trna")]
t1=sampleSubset[str_detect(sampleSubset,"tp.1")]
t4=sampleSubset[str_detect(sampleSubset,"tp.4")]
sampleIDsRNA=c(t1,t4)

dirsSubset=kallistoDirs[str_detect(kallistoDirs,"trna")]
t1=dirsSubset[str_detect(dirsSubset,"tp.1")]
t4=dirsSubset[str_detect(dirsSubset,"tp.4")]
kallistoDirsRNA=c(t1,t4)

# 2. handle metadata reading and path association
s2cRBF=read.table(metadataFileRBF,header=TRUE,stringsAsFactors=FALSE,sep=",")
s2cRBF=dplyr::mutate(s2cRBF,path=kallistoDirsRBF)

s2cRNA=read.table(metadataFileRNA,header=TRUE,stringsAsFactors=FALSE,sep=",")
s2cRNA=dplyr::mutate(s2cRNA,path=kallistoDirsRNA)

# 3. running Sleuth
# 3.1. running Sleuth for RBF
soRBF=sleuth_prep(s2cRBF,extra_bootstrap_summary=TRUE,num_cores=1)
soRBF=sleuth_fit(soRBF,~condition,'full')
soRBF=sleuth_fit(soRBF,~1,'reduced')
soRBF=sleuth_lrt(soRBF,'reduced','full')

sleuthTable=sleuth_results(soRBF,'reduced:full','lrt',show_all=FALSE)
sleuthSignificant=dplyr::filter(sleuthTable,qval<=0.05)
dim(sleuthSignificant)
#head(sleuthSignificant,20)
#plot_bootstrap(soRBF, "lcl|AE004437.1_gene_284", units = "est_counts", color_by = "condition")
write.table(sleuthSignificant,file=sleuthResultsRBF,sep=",")

# 3.2. running Sleuth for RNA
soRNA=sleuth_prep(s2cRNA,extra_bootstrap_summary=TRUE,num_cores=1)
soRNA=sleuth_fit(soRNA,~condition,'full')
soRNA=sleuth_fit(soRNA,~1,'reduced')
soRNA=sleuth_lrt(soRNA,'reduced','full')

sleuthTable=sleuth_results(soRNA,'reduced:full','lrt',show_all=FALSE)
sleuthSignificant=dplyr::filter(sleuthTable,qval<=0.05)
dim(sleuthSignificant)
#head(sleuthSignificant,20)
#plot_bootstrap(soRNA, "lcl|AE004437.1_gene_1082", units = "est_counts", color_by = "condition")
write.table(sleuthSignificant,file=sleuthResultsRNA,sep=",")
