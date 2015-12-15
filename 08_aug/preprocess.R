#!/usr/bin/env Rscript

library(dplyr)

min.otu.size <- 1000
max.blank.frac <- 0.10
max.sediment.frac <- 0.05

# read in OTU table
otu <- read.table("otu_raw.txt", header=T, sep="\t", row.names=1)

# pool replicate samples, the drop the originals
otu$M3 <- otu$M3.1 + otu$M3.2
otu$M8 <- otu$M8.1 + otu$M8.2
otu <- otu[, !(colnames(otu) %in% c("M3.1", "M3.2", "M8.1", "M8.2"))]

# note which OTUs that have fewer than 1000 counts
too.small.otus <- rowSums(otu) < min.otu.size

# normalize by column (i.e., convert to relative abundances)
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

# remove OTUs with small counts
otu <- otu[!too.small.otus, ]

# remove OTUs w/ more than 10% of their reads in the blanks
otu <- otu[otu$MEB + otu$MSB <= max.blank.frac, ]

# remove OTUs w/ more than 5% of their reads in the sediment-y sample (M22)
otu <- otu[otu$M22 <= max.sediment.frac, ]

# remove extra columns
otu <- otu[, !(colnames(otu) %in% c("MEB", "MSB", "M22"))]

write.table(otu, "otu.txt", sep="\t", row.names=T, quote=F)


## postprocess
# prepare a data frame for plotting
otu$id <- names(groups)
otu$group <- groups
otu <- melt(otu, id.vars=c("id", "group"), variable.name="sample")
otu$depth <- sapply(otu$sample %>% as.list, function(x) sub("M", "", x) %>% as.numeric)
write.table(otu, file=plot.dat.fn, sep="\t", quote=F)
