### Figure 5A
CorPlot <- function(genelist, table = sepOnlyMin15counts) {
  require(ggplot2)
  data <- table[table$gene %in% genelist,]
  d1 <- data[data$timepoint == '1',]
  d1$panel <- 'time point 1'
  d2 <- data[data$timepoint == '2',]
  d2$panel <- 'time point 2'
  d3 <- data[data$timepoint == '3',]
  d3$panel <- 'time point 3'
  d4 <- data[data$timepoint == '4',]
  d4$panel <- 'time point 4'
  data <- rbind(d1,d2,d3,d4)
  
  ggplot(data, aes(x = mrnaCount, y = rpfCount)) +
    facet_grid(panel~.,scales = 'fixed') + 
    theme(text = element_text(size=20)) + 
    geom_point(alpha = .6, size = 2.5) +
    geom_smooth(method = 'lm', se = FALSE) +
    scale_x_log10(limits = c(10,10000000)) +
    scale_y_log10(limits = c(10,100000)) +
    ylab("RPF normalized count") +
    xlab("mRNA normalized count")
}

multiplot(CorPlot(unique(sepOnlyMin15counts$gene)), CorPlot(riboVNGm), cols = 2)

### Figure 5B

teBoxplot <- function(genelist) {
  require(ggplot2)
  data <- meltedCombined[meltedCombined$gene %in% genelist & meltedCombined$sampleType %in% 'TE',]
    
  ggplot(data, aes(x = gene, y = count, group = interaction(timepoint,gene))) +
    theme(text = element_text(size=20)) + 
    geom_boxplot(aes(color = timepoint)) +
    ylab("RPF/mRNA count ratio") +
    xlab("gene") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_log10(limits=c(.10,10)) +
    scale_x_discrete(labels = geneNames[genelist,2])
}

teBoxplot(c('VNG1104G','VNG1105G','VNG1137G','VNG1139Gm','VNG1691G','VNG1714G','VNG2514G'))