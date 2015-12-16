#!/usr/bin/env Rscript

library(dplyr)

otu <- read.table('otu.txt', header=T, row.names=1, sep="\t")
depths <- sapply(colnames(otu), function(x) sub("M", "", x) %>% as.numeric)
mean.depth <- apply(otu, 1, function(x) weighted.mean(depths, w=x))
dat <- data.frame(otu=rownames(otu), mean.depth=mean.depth)
write.table(dat, "mean-depths.txt", sep="\t", row.names=F, quote=F)
