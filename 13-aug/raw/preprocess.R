#!/usr/bin/env Rscript

library(dplyr)

min.otu.size <- 100
max.sediment.frac <- 0.05

# read in OTU table
otu <- read.table("otu_rename.txt", header=T, sep="\t", row.names=1)

sum.columns <- function(otu, res, cols) {
  otu[[res]] <- 0
  for (col in cols) otu[[res]] <- otu[[res]] + otu[[col]]
  otu <- otu[, !(colnames(otu) %in% cols)]
  otu
}

# pool replicate samples, the drop the originals
otu <- otu %>%
  sum.columns("M10", c("M10.1", "M10.2")) %>%
  sum.columns("M11", c("M11.1", "M11.2")) %>%
  sum.columns("M14", c("M14.1", "M14.2")) %>%
  sum.columns("M17", c("M17.1", "M17.2")) %>%
  sum.columns("M2", c("M2.1", "M2.2")) %>%
  sum.columns("M7", c("M7.1", "M7.2")) %>%
  sum.columns("M9", c("M9.1", "M9.2"))

too.small.otus <- rowSums(otu) < min.otu.size

# normalize by column (i.e., convert to relative abundances)
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

# remove OTUs with small counts
otu <- otu[!too.small.otus, ]

# remove OTUs w/ more than 5% of their reads in the sediment-y sample (M22)
otu <- otu[otu$M22 <= max.sediment.frac, ]

# remove extra columns
otu <- otu[, !(colnames(otu) %in% c("M22"))]

# put the OTU IDs back
sample.ids <- colnames(otu)
otu$OTU_ID <- rownames(otu)
otu <- otu[, c('OTU_ID', sample.ids)]

write.table(otu, "otu.txt", sep="\t", row.names=T, quote=F)
