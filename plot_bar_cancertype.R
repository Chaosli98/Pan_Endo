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

####color settings
ecol = list()
ecol[['cellType.major']]=c("#7FC97F","#386CB0","#FDC086","#BF5B17","#BEAED4","#666666")
names(ecol[['cellType.major']])=c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics')


meta <- read.table('../data/panE.freq.cancerTissue.majorC.csv',
                 header = TRUE,stringsAsFactors = FALSE,sep = ',')
meta <- separate(data = meta, col = X, into = c("cancerType", "tissue"), sep = "\\|")
meta <- melt(meta,id=c('cancerType','tissue'))


meta$tissue = factor(meta$tissue,levels = c('NT','T'))
meta$variable = factor(meta$variable,levels = c('tip.cell','arteries','veins','capillaries','hypoxia','lymphatics'))
meta$value <- as.numeric(meta$value )

##fig 5B

p <- ggplot(meta, aes(x=tissue, y=value, fill=variable)) +
  geom_bar(stat="identity", position="fill", colour="white", linewidth=0.5, width=0.8) +
  scale_fill_manual(values=ecol$cellType.major) + 
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  facet_wrap(facets = 'cancerType',nrow = 2)
p
ggsave('../out/major/eachcancerType.tissue.major.bar.pdf',plot = p,width = 9,height = 4)













