#!/usr/bin/env Rscript

library(dplyr)

min.otu.size <- 1000
max.blank.frac <- 0.10
max.sediment.frac <- 0.05

# read in OTU table
otu <- read.table("raw_rename.txt", header=T, sep="\t", row.names=1)

sum.columns <- function(otu, res, cols) {
  otu[[res]] <- 0
  for (col in cols) otu[[res]] <- otu[[res]] + otu[[col]]
  otu <- otu[, !(colnames(otu) %in% cols)]
  otu
}

# pool replicate samples, the drop the originals
otu <- sum.columns(otu, "M17", c("M17.1", "M17.2", "M17.3"))
otu <- sum.columns(otu, "M03", c("M03.1", "M03.2", "M03.4", "M03.5", "M03.6", "M03.7", "M03.8"))
otu <- sum.columns(otu, "M08", c("M08.1", "M08.2"))
otu <- sum.columns(otu, "M09", c("M09.1", "M09.2"))

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
