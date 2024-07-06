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

###only use T-NT paired samples
df = read.table('../out_TNpaired/Counter.TNp.csv',
                header = TRUE,stringsAsFactors = FALSE,sep = ',')
rownames(df) = df$cellType.sub
df = df[-1]
ch = chisq.test(df)

ch$residuals

heatdf = data.frame(ch$residuals)
rownames(heatdf) = rownames(df)


ecol = readRDS('../data/panE.color.rds')

df2 = heatdf
df2$subset = rownames(heatdf)
df2$major = c('tip.cell','tip.cell','tip.cell','tip.cell','tip.cell',
              'veins','veins','veins','veins','veins',
              'capillaries','capillaries','capillaries',
              'arteries','arteries','arteries',
              'hypoxia','lymphatics')


pdf('../out_TNpaired/Counter.TNp.lolipop.sorting.pdf',width = 8,height = 5)
ggdotchart(df2, y = "T", x = "subset",
           color = "major",    
           sorting = "descending",                        
           add = "segments",                             
           xlab="", 
           rotate = FALSE,
           #group = "type", 
           dot.size = 5,
           palette = ecol$cellType.major
)
dev.off()
