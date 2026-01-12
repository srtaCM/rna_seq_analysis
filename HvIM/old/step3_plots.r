root <- 0
sig2 <- 1
mu <- 0
setwd("~/Desktop/GenFinal")
SixSamples.readcount <- read.delim("~/Desktop/GenFinal/SixSamples.readcount.out", comment.char="#")
#View(SixSamples.readcount)
library(tidyverse)
library(dplyr)
Samples<- SixSamples.readcount %>% select(1, 6:12)
View(Samples)
install.packages("matrixStats")
library(matrixStats)
Samples$row_std = rowSds(as.matrix(Samples[c(3:8)]))
View(Samples)
samples2<- Samples %>% arrange(desc(row_std)) %>% slice_head(n = 2000)
View(samples2)
library("ggplot2")
heatsample<- data.matrix(samples2 %>% select(3:8))
heatsample2<- log2(heatsample+1)
summary(heatsample)
View(heatsample2)
heatmap(heatsample2, Colv = NULL, Rowv = NULL, scale="row", col = cm.colors(5000), xlab="Length", ylab="rc[0g2", main="heatmap")


scattersample<- data.frame(samples2 %>% select(2:8))
View(scattersample)
scattersample2<- log2(scattersample+1)
View(scattersample2)
ggplot(scattersample2, aes(x=Length , y=SRR2927328Aligned.sortedByCoord.out.bam)) + geom_point(size=2, shape=23)

