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

#######color
ecol = list()
ecol[['cellType.major']]=c("#7FC97F","#386CB0","#FDC086","#BF5B17","#BEAED4","#666666")
names(ecol[['cellType.major']])=c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics')
ecol[['tissue']]=c('#1F78B4','#FE9E37')
names(ecol[['tissue']])=c('NT','T')

#######read data
meta <- read.table('../data/panE.freq.eachsample.majorC.csv',
                   header = TRUE,stringsAsFactors = FALSE,sep = ',')
sample <- read.table('../data/panE.sample.csv',
                     header = TRUE,stringsAsFactors = FALSE,sep = ',')
rownames(meta) = meta$X
meta = meta[sample$sampleID,]


meta$tissue = sample$tissue
meta$cancerType = sample$cancerType
meta <- melt(meta,id=c('X','tissue','cancerType'))

meta$tissue = factor(meta$tissue,levels = c('NT','T'))
meta$variable = factor(meta$variable,levels = c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics'))
meta$value <- as.numeric(meta$value )

##fig 5A


order_vec <- setNames(c('NT','T'),c('NT','T'))
pp.list <- llply(order_vec, function(i){
  df = meta[meta$tissue==i,]
  p <- ggplot(df, mapping = aes(
    x = variable, 
    y = value,
    color = variable)) + 
    geom_boxplot(width=0.7, alpha=0.5, size=0.7,outlier.shape = NA) + 
    geom_quasirandom(size=0.5, alpha=0.7, shape=16, width=0.3)+
    theme_classic() +
    theme(legend.position="none", panel.background = element_rect(fill = 'white',color = "black")) + 
    stat_compare_means() + 
    scale_color_manual(values=ecol$cellType.major) +
    ggtitle(i)
}, .parallel=T)
pp <- plot_grid(plotlist=pp.list,align="hv",ncol=2)
pp
ggsave(pp, file='../out/major/eachsample.tissue.major.box.pdf', width=7, height=2.8, limitsize=F)



order_vec <- setNames(c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics'),
                      c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics'))
pp.list <- llply(order_vec, function(i){
  df = meta[meta$variable==i,]
  p <- ggplot(df, mapping = aes(
    x = tissue, 
    y = value,
    color = tissue)) + 
    geom_boxplot(width=0.7, alpha=0.5, size=0.7,outlier.shape = NA) + 
    geom_quasirandom(size=0.5, alpha=0.7, shape=16, width=0.3)+
    theme_classic() +
    theme(legend.position="none", panel.background = element_rect(fill = 'white',color = "black")) + 
    stat_compare_means(method = "wilcox.test",size=2.5) + 
    scale_color_manual(values=c('#1F78B4','#FE9E37')) +
    ggtitle(i)
}, .parallel=T)
pp <- plot_grid(plotlist=pp.list,align="hv",ncol=6)
pp

ggsave(pp, file='../out/major/eachsample.major.tissue.box.pdf', width=11, height=2.8, limitsize=F)

