library(tidyverse)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggsci)
library(cowplot)

#install.packages("vegan")
library(vegan)

ecol = readRDS('../data/panE.color.rds')

############# fig 4D
meta = read.table('../data/panE.freq.eachcancerType.subC.csv',
                  header = TRUE,stringsAsFactors = FALSE,sep = ',')
meta <- melt(meta,id=c('X'))
meta$variable = factor(meta$variable,levels = c('E01.tip.KDR','E02.tip.CXCR4','E03.tip.DNAJB1','E04.tip.ATP5E',
                                                'E05.tip.STMN1','E06.veins.SELE','E07.veins.CLU','E09.veins.COL4A1',
                                                'E08.veins.FABP4','E10.veins.IGFBP5','E11.capillaries.TMEM100',
                                                'E12.capillaries.FABP4','E13.capillaries.RGCC','E14.arteries.DNAJB1',
                                                'E15.arteries.CLU','E16.arteries.COL1A1','E17.hypoxia.MT1X',
                                                'E18.lymphatics.PROX1'))

meta = meta[meta$X!='TGCT',]
meta = meta[meta$X!='CESC',]

meta2 <- dcast(meta,X~variable)
rownames(meta2) <- meta2$X
meta2 <- meta2[,-1]
meta2[is.na(meta2)] <- 0


tree = hclust(vegan::vegdist(meta2, method = 'bray'), 
              method = 'average')

p1 = ggtree(tree) + geom_tiplab() +xlim(NA,2)
p1


g2 <- ggplot(meta,aes(y=value,x=X,fill=variable,colour = X))+
  geom_bar(stat = "identity",width = 0.5,alpha=1,color=NA) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank())+
  
  scale_x_discrete(limits = rev(c('STAD','CRC','NPC','PAAD','ESCA','PRAD','LUAD','BRCA','THCA','CHOL','SKCM',
                                  'LIHC','OV','GM','RCC','HNSC','cSCC')),) +
  
  scale_fill_manual(values = ecol$cellType.sub) +
  scale_colour_manual(values = ecol$cellType.sub) +
  
  coord_flip()+
  ylab("Frequency") +
  xlab("") 
g2


g <- ggdraw()+
  draw_plot(p1, 0, 0.05, 0.7, 0.95)+
  draw_plot(g2, 0.2, 0, 0.8, 1)
g

ggsave("../out/tree/sub.pdf", g, height = 5, width =9)

