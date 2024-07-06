suppressPackageStartupMessages({
  library("plyr")
  library("dplyr")
  library("data.table")
  library("reshape2")
  library("tibble")
  library("ggsci")
  library("ggpubr")
  library("ggrastr")
  library("ggbeeswarm")
  library("cowplot")
  library("stringr")
  library("future")
  library("purrr")
  library("furrr")
  library("ComplexHeatmap")
  library("circlize")
})

library(tidyverse)
library(ggplot2)

###related to fig S1F, 4E and S4E

ecol = readRDS(file = '../data/panE.color.rds')


meta <- read.table('../data/panE.freq.eachsample.subC.csv',
                   header = TRUE,stringsAsFactors = FALSE,sep = ',')
sample <- read.table('../data/panE.sample.csv',
                     header = TRUE,stringsAsFactors = FALSE,sep = ',')

rownames(meta) = meta$X
meta = meta[sample$sampleID,]
meta$tissue = sample$tissue
meta$cancerType = sample$cancerType
meta = meta[meta$tissue=='T',]
meta <- melt(meta,id=c('X','tissue','cancerType'))

meta$variable = factor(meta$variable,levels = c('E01.tip.KDR','E02.tip.CXCR4','E03.tip.DNAJB1','E04.tip.ATP5E',
                                                'E05.tip.STMN1','E06.veins.SELE','E07.veins.CLU','E09.veins.COL4A1',
                                                'E08.veins.FABP4','E10.veins.IGFBP5','E11.capillaries.TMEM100',
                                                'E12.capillaries.FABP4','E13.capillaries.RGCC','E14.arteries.DNAJB1',
                                                'E15.arteries.CLU','E16.arteries.COL1A1','E17.hypoxia.MT1X',
                                                'E18.lymphatics.PROX1'))
meta$value <- as.numeric(meta$value )
meta$cancerType <- factor(meta$cancerType,levels = sort(unique(meta$cancerType)))

order_vec <- setNames(c('E01.tip.KDR','E02.tip.CXCR4','E03.tip.DNAJB1','E04.tip.ATP5E',
                        'E05.tip.STMN1','E06.veins.SELE','E07.veins.CLU','E09.veins.COL4A1',
                        'E08.veins.FABP4','E10.veins.IGFBP5','E11.capillaries.TMEM100',
                        'E12.capillaries.FABP4','E13.capillaries.RGCC','E14.arteries.DNAJB1',
                        'E15.arteries.CLU','E16.arteries.COL1A1','E17.hypoxia.MT1X',
                        'E18.lymphatics.PROX1'),
                      c('E01.tip.KDR','E02.tip.CXCR4','E03.tip.DNAJB1','E04.tip.ATP5E',
                        'E05.tip.STMN1','E06.veins.SELE','E07.veins.CLU','E09.veins.COL4A1',
                        'E08.veins.FABP4','E10.veins.IGFBP5','E11.capillaries.TMEM100',
                        'E12.capillaries.FABP4','E13.capillaries.RGCC','E14.arteries.DNAJB1',
                        'E15.arteries.CLU','E16.arteries.COL1A1','E17.hypoxia.MT1X',
                        'E18.lymphatics.PROX1'))

meta = meta[meta$cancerType!='TGCT',]
meta = meta[meta$cancerType!='CESC',]

pp.list <- llply(order_vec, function(i){
  df = meta[meta$variable==i,]
  med = aggregate(df$value, by=list(cancerT=df$cancerType),median)
  df$cancerType <- factor(df$cancerType,levels = rev(med$cancerT[order(med$x)]))
  p <- ggplot(df, mapping = aes(
    x = cancerType, 
    y = value,
    color = cancerType)) + 
    geom_boxplot(width=0.7, alpha=0.5, size=0.5,outlier.shape = NA) + 
    geom_quasirandom(size=0.5, alpha=0.7, shape=16, width=0.3)+
    theme_classic() +
    theme(legend.position="none", panel.background = element_rect(fill = 'white',color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1)) + 
    stat_compare_means(size=4) + 
    scale_color_manual(values=ecol$cancerType) +
    ggtitle(i)+
    xlab('')+
    ylab('Proportion')
}, .parallel=T)
pp <- plot_grid(plotlist=pp.list,align="hv",ncol=2)
pp
ggsave(pp, file='../out/sub/eachsub.cancerType.box.pdf', width=12, height=20, limitsize=F)



















