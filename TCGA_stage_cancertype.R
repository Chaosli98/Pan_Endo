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
library(Rmisc)

adf = read.table('../data/TCGA/TCGA_ssgsea.score.csv', header=T, stringsAsFactors=F, check.names=F, sep=",")

dat = read.table("../data/TCGA/TCGA_surdf.csv", header=T, stringsAsFactors=F, check.names=F, sep=",")

rownames(dat)=dat$sample
rownames(adf)=adf$sample
use=intersect(rownames(dat),rownames(adf))
dat = dat[use,]
adf = adf[use,]


dat = cbind(dat,adf)
dat$cancerType = dat$cancer.type.abbreviation
dat = dat[!is.na(dat$ajcc_pathologic_tumor_stage),]
dat$stage = dat$ajcc_pathologic_tumor_stage
unique(dat$stage)

dat = dat[grepl("Stage I",dat$stage),]
unique(dat$stage)

dat$stage = gsub("[ABCD]$","",dat$stage)
unique(dat$stage)



ecol = readRDS('../data/panE.color.rds')
ecol$cellType.sub['E02.tip.CXCR4']
ecol$cellType.sub['E06.veins.SELE']


new = c()
unique(dat$stage)

for (each in dat$stage) {
  if (each %in% c('Stage I','Stage II')) {
    new = append(new,'early')
  }else{
    new = append(new,'late')
  }
}
dat$phase = new



##########SEM
pp.list <- llply(order_vec, function(i){
  ddat = dat[dat$cancer.type.abbreviation==i,]
  
  ddat2 = ddat[,c('phase','E06.veins.SELE','E02.tip.CXCR4')]
  ddat2 = melt(ddat2,id=c('phase'))
  colnames(ddat2) = c('stage','type','score')
  tgc <- summarySE(ddat2, measurevar="score", groupvars=c("type","stage"))
  
  
  df1 = aggregate(ddat$E06.veins.SELE, by=list(stage=ddat$phase),mean)
  df1$type = 'E06.veins.SELE'
  df3 = aggregate(ddat$E02.tip.CXCR4, by=list(stage=ddat$phase),mean)
  df3$type = 'E02.tip.CXCR4'
  
  df = rbind(df1,df3)
  df$se = tgc$se
  df$stage = factor(df$stage,levels = c('early','late'))
  df$type = factor(df$type,levels = c('E06.veins.SELE','E02.tip.CXCR4'))
  
  g = ggplot(data = df, mapping = aes(x = stage, y = x, group = type, colour = type, fill = type,shape=type))+
    geom_line(size=1,aes(linetype=type)) + 
    geom_point(size=1.5) + 
    scale_linetype_manual(values = c('solid','twodash'))+ 
    scale_color_manual(values = c("#89CB72","#FD9D37"))+
    scale_shape_manual(values = c(21,23))+
    scale_fill_manual(values = c("#89CB72","#FD9D37"))+
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
    xlab('Tumor stage')+
    ylab('score')+
    ggtitle(i)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_grid(type~.,scales = 'free')+
    geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.2)
}, .parallel=T)

pp <- plot_grid(plotlist=pp.list,align="hv",ncol=3)

ggsave(filename = '../out/TCGA_SELE.CXCR4.score.stage.percancer.facet.sem.phase.pdf',pp,width = 11,height = 28)
