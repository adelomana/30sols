###Figure 6A

CorPlot <- function(genelist, replicate = 1:3, table = tmp3) {
  require(ggplot2)
  data <- table[table$gene %in% genelist & table$replicate %in% replicate,]
  
  rpfData <- data[data$sampleType == 'RPF',1:4]
  rpfData <- rbind(rpfData, tmp5[tmp5$sampleType == 'RPF', 1:4])
  rpfData$gene <- gsub('rpf','avg',rpfData$gene)
  names(rpfData)[2] <- 'RPFlog2_sfc'
  
  riboData <- data[data$sampleType == 'ribosome',1:4]
  riboData <- rbind(riboData, tmp5[tmp5$sampleType == 'ribosome', 1:4])
  riboData$gene <- gsub('ribosome','avg',riboData$gene)
  names(riboData)[2] <- 'Ribolog2_sfc'
  
  
  data <- merge(rpfData,riboData)
  data$category <- annotations[data$gene,1]
  
  d2 <- data[data$timepoint == '2',]
  d2$panel <- 'time point 2'
  d3 <- data[data$timepoint == '3',]
  d3$panel <- 'time point 3'
  d4 <- data[data$timepoint == '4',]
  d4$panel <- 'time point 4'
  data <- rbind(d2,d3,d4)
   
  ggplot(data, aes(RPFlog2_sfc, Ribolog2_sfc, group = interaction(replicate,gene))) +
    theme(text = element_text(size=20)) + 
    facet_grid(panel~.,scales = 'fixed') + 
    geom_point(aes(color = category, shape = category), alpha = .35, size = 3) +
    #geom_path(aes(color = gene), arrow = arrow(angle = 15, type = 'closed')) + 
    ggtitle('correlation over time') +
    #coord_flip() + 
    xlim(-7, 3) +
    ylim(-4, 2) +
    xlab("RPF log2 fold change") +
    ylab("riboproteome log2 fold change")
}

CorPlot2 <- function(genelist, replicate = 1:3, table = tmp3) {
  require(ggplot2)
  data <- table[table$gene %in% genelist & table$replicate %in% replicate,]
  
  mrnaData <- data[data$sampleType == 'mRNA',1:4]
  mrnaData <- rbind(mrnaData, tmp5[tmp5$sampleType == 'mRNA', 1:4])
  mrnaData$gene <- gsub('mrna','avg',mrnaData$gene)
  names(mrnaData)[2] <- 'mrnalog2_sfc'
  
  lysateData <- data[data$sampleType == 'lysate',1:4]
  lysateData <- rbind(lysateData, tmp5[tmp5$sampleType == 'lysate', 1:4])
  lysateData$gene <- gsub('lysate','avg',lysateData$gene)
  names(lysateData)[2] <- 'lysatelog2_sfc'
  
  data <- merge(mrnaData,lysateData)
  data$category <- annotations[data$gene,1]
  
  d2 <- data[data$timepoint == '2',]
  d2$panel <- 'time point 2'
  d3 <- data[data$timepoint == '3',]
  d3$panel <- 'time point 3'
  d4 <- data[data$timepoint == '4',]
  d4$panel <- 'time point 4'
  data <- rbind(d2,d3,d4)
  
  ggplot(data, aes(mrnalog2_sfc, lysatelog2_sfc, group = interaction(replicate,gene))) +
    theme(text = element_text(size=20)) + 
    facet_grid(panel~.,scales = 'fixed') + 
    geom_point(aes(color = category, shape = category), alpha = .35, size = 3) +
    #geom_path(aes(color = gene), arrow = arrow(angle = 15, type = 'closed')) + 
    ggtitle('correlation over time') +
    #coord_flip() + 
    xlim(-7, 3) +
    ylim(-4, 2) +
    xlab("mRNA log2 fold change") +
    ylab("lysate proteome log2 fold change")
}

multiplot(CorPlot(c(riboVNG,eifVNG,aarsVNG,rnaPolVNG)), CorPlot2(c(riboVNG,eifVNG,aarsVNG,rnaPolVNG)), cols = 2)

### Figure 6B

testCorPlot2_mrna <- function(genelist, replicate = 1:3, table = tmp3) {
  require(ggplot2)
  data <- table[table$gene %in% genelist & table$replicate %in% replicate,]
  
  mrnaData <- data[data$sampleType == 'mRNA',1:3]
  mrnaFilData <- data.frame(gene = character(), log2_sfc = numeric(), timepoint = integer())
  for(g in genelist) {
    for (tp in 1:4) {
      temp <- mean(mrnaData[mrnaData$gene %in% g & mrnaData$timepoint == tp, 'log2_sfc'])
      mrnaFilData <- rbind(mrnaFilData, data.frame(gene = g, log2_sfc = temp, timepoint = tp))
    }
  }
  mrnaFilData <- rbind(mrnaFilData, tmp5[tmp5$sampleType == 'mRNA', 1:3])
  mrnaFilData$gene <- gsub('mrna','RP median',mrnaFilData$gene)
  names(mrnaFilData)[2] <- 'mrnalog2_sfc'
  
  lysateData <- data[data$sampleType == 'lysate',1:3]
  lysateFilData <- data.frame(gene = character(), log2_sfc = numeric(), timepoint = integer())
  for(g in genelist) {
    for (tp in 1:4) {
      temp <- mean(lysateData[lysateData$gene %in% g & lysateData$timepoint == tp, 'log2_sfc'])
      lysateFilData <- rbind(lysateFilData, data.frame(gene = g, log2_sfc = temp, timepoint = tp))
    }
  }
  lysateFilData <- rbind(lysateFilData, tmp5[tmp5$sampleType == 'lysate', 1:3])
  lysateFilData$gene <- gsub('lysate','RP median',lysateFilData$gene)
  names(lysateFilData)[2] <- 'lysatelog2_sfc'
  
  
  data <- merge(mrnaFilData,lysateFilData)
  data$timepoint <- factor(data$timepoint)
  
  ggplot(data, aes(mrnalog2_sfc, lysatelog2_sfc, group = gene)) +
    theme(text = element_text(size=20)) + 
    geom_point(aes(color = gene, shape = timepoint), alpha = .35, size = 3) +
    geom_path(aes(color = gene), arrow = arrow(angle = 15, type = 'closed'), alpha = .35) + 
    xlim(-5, 1) +
    ylim(-2, 2) +
    xlab("mRNA log2 fold change") +
    ylab("lysate log2 fold change")
}

testCorPlot2_rpf <- function(genelist, replicate = 1:3, table = tmp3) {
  require(ggplot2)
  data <- table[table$gene %in% genelist & table$replicate %in% replicate,]
  
  mrnaData <- data[data$sampleType == 'RPF',1:3]
  mrnaFilData <- data.frame(gene = character(), log2_sfc = numeric(), timepoint = integer())
  for(g in genelist) {
    for (tp in 1:4) {
      temp <- mean(mrnaData[mrnaData$gene %in% g & mrnaData$timepoint == tp, 'log2_sfc'])
      mrnaFilData <- rbind(mrnaFilData, data.frame(gene = g, log2_sfc = temp, timepoint = tp))
    }
  }
  mrnaFilData <- rbind(mrnaFilData, tmp5[tmp5$sampleType == 'RPF', 1:3])
  mrnaFilData$gene <- gsub('rpf','RP median',mrnaFilData$gene)
  names(mrnaFilData)[2] <- 'mrnalog2_sfc'
  
  lysateData <- data[data$sampleType == 'ribosome',1:3]
  lysateFilData <- data.frame(gene = character(), log2_sfc = numeric(), timepoint = integer())
  for(g in genelist) {
    for (tp in 1:4) {
      temp <- mean(lysateData[lysateData$gene %in% g & lysateData$timepoint == tp, 'log2_sfc'])
      lysateFilData <- rbind(lysateFilData, data.frame(gene = g, log2_sfc = temp, timepoint = tp))
    }
  }
  lysateFilData <- rbind(lysateFilData, tmp5[tmp5$sampleType == 'ribosome', 1:3])
  lysateFilData$gene <- gsub('ribosome','RP median',lysateFilData$gene)
  names(lysateFilData)[2] <- 'lysatelog2_sfc'
  
  
  data <- merge(mrnaFilData,lysateFilData)
  data$timepoint <- factor(data$timepoint)
  
  ggplot(data, aes(mrnalog2_sfc, lysatelog2_sfc, group = gene)) +
    theme(text = element_text(size=20)) + 
    geom_point(aes(color = gene, shape = timepoint), alpha = .35, size = 3) +
    geom_path(aes(color = gene), arrow = arrow(angle = 15, type = 'closed'), alpha = .35) + 
    xlim(-5, 1) +
    ylim(-2, 2) +
    xlab("RPF log2 fold change") +
    ylab("riboproteome log2 fold change")
}


multiplot(testCorPlot2_rpf(c('VNG0433C','VNG1158G')), testCorPlot2_mrna(c('VNG0433C', 'VNG1158G')), cols = 2)