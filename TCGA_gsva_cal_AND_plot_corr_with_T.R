library(AUCell)
library(Seurat)
library(SingleCellExperiment)
suppressPackageStartupMessages({
  library("plyr")
  library("dplyr")
  library("ggpubr")
  library("ggsci")
  library("reshape2")
  library("survival")
  library("survminer")
  library("data.table")
})
library(GSVA)
library(ggbeeswarm)


sce = readRDS('../data/TCGA/TCGA_Toil.rds')
sce = sce[,sce$sampleType=="Primary.Solid.Tumor"]
assay(sce,'exprs') = log1p(assay(sce,'tpm'))
se = as.Seurat(sce,counts=NULL,data='exprs')

####cal score for each cancer type
geneSets = list('T'=c('CD3D','CD3E','CD3G'),
                'CD8T'=c('CD8A','CD8B'),
                'E02.tip.CXCR4'=c('IGFBP3', 'ANGPT2', 'INSR', 'ITGB1', 'UNC5B', 'ESM1', 'CXCR4', 'SPARC'),
                'E06.veins.SELE'=c('SELE', 'ICAM1', 'CCL2', 'RND1', 'NFKBIA', 'EIF1', 'SOD2', 'CXCL2')
                )


new = 'True'
for (eachc in unique(se$cancer.type.abbreviation)) {
  se.c = subset(se,cancer.type.abbreviation==eachc)
  
  tmp = gsva(as.matrix(se.c@assays$originalexp@data), geneSets, method="ssgsea", kcdf="Gaussian", ssgsea.norm=T, parallel.sz=40)
  
  
  bdf = data.frame(t(tmp))
  bdf$sample=se.c$sample
  bdf$cancer.type.abbreviation=eachc

  if (new=='True') {
    adf = bdf
    new='added'
  }else{
    adf = rbind(adf,bdf)
  }
}


adf
write.table(adf, "../data/TCGA/TCGA_ssgsea.score.csv", sep=",", quote=F, col.names=T, row.names=F)

##################################################corr  fig 2G and S3D
adf = read.table('../data/TCGA/TCGA_ssgsea.score.csv', header=T, stringsAsFactors=F, check.names=F, sep=",")

colnames(adf)

order_vec <- setNames(unique(adf$cancer.type.abbreviation),unique(adf$cancer.type.abbreviation))

pp.list <- llply(order_vec, function(i){
  df = adf[adf$cancer.type.abbreviation==i,]
  
  p <- ggplot(df, aes(x=E06.veins.SELE, y=T)) +
    geom_point(shape=16,size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm',formula=y~x) +
    theme_classic2() +
    theme(strip.text.x=element_text(size=8), legend.position="none")+
    ggtitle(i)
}, .parallel=T)
pp <- plot_grid(plotlist=pp.list,align="hv",ncol=4)
ggsave(filename = '../out/TCGA_ssgsea.cancertype.SELE.T.pdf',pp,width = 10,height = 20)



pp.list <- llply(order_vec, function(i){
  df = adf[adf$cancer.type.abbreviation==i,]
  
  p <- ggplot(df, aes(x=E06.veins.SELE, y=CD8T)) +
    geom_point(shape=16,size=1) +
    stat_cor(size=3) +
    geom_smooth(method='lm',formula=y~x) +
    theme_classic2() +
    theme(strip.text.x=element_text(size=8), legend.position="none")+
    ggtitle(i)
}, .parallel=T)
pp <- plot_grid(plotlist=pp.list,align="hv",ncol=4)
ggsave(filename = '../out/TCGA_ssgsea.cancertype.SELE.CD8T.pdf',pp,width = 10,height = 20)










