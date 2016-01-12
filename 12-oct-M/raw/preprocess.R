#!/usr/bin/env Rscript

library(dplyr)

min.otu.size <- 0.001
max.blank.frac <- 0.10
max.sediment.frac <- 0.05

# read in OTU table
otu <- read.table("otu_counts.txt", header=T, sep="\t", row.names=1)

sum.columns <- function(otu, res, cols) {
  otu[[res]] <- 0
  for (col in cols) otu[[res]] <- otu[[res]] + otu[[col]]
  otu <- otu[, !(colnames(otu) %in% cols)]
  otu
}

rename.columns <- function(otu, res, cols) {
  otu[[res]] <- 0
  otu[[res]] <- otu[[cols]]
  otu <- otu[, !(colnames(otu) %in% cols)]
  otu
}

# pool replicate samples, the drop the originals
otu <- rename.columns(otu, "M01", c("SB100212TAWMD01VV4TMR1"))
otu <- rename.columns(otu, "M03", c("SB100212TAWMD03VV4TMR1"))
otu <- rename.columns(otu, "M05", c("SB100212TAWMD05VV4TMR1"))
otu <- rename.columns(otu, "M07", c("SB100212TAWMD07VV4TMR2"))
otu <- rename.columns(otu, "M09", c("SB100212TAWMD09VV4TMR1"))
otu <- rename.columns(otu, "M11", c("SB100212TAWMD11VV4TMR2"))
otu <- rename.columns(otu, "M15", c("SB100212TAWMD15VV4TMR2"))
otu <- rename.columns(otu, "M17", c("SSB100212TAWMD17VV4TMR2"))
otu <- rename.columns(otu, "M19", c("SB100212TAWMD19VV4TMR1"))

# normalize by column (i.e., convert to relative abundances)
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

# remove OTUs w/ more than 10% of their reads in the blanks
#otu <- otu[otu$MEB  <= max.blank.frac, ]

# remove OTUs w/ more than 5% of their reads in the sediment-y sample (M22)
# otu <- otu[otu$M22 <= max.sediment.frac, ]

# remove extra columns
otu <- otu[, !(colnames(otu) %in% c("SB100212TAWMD06VV4TMR1", "SB100212TAWMD08VV4TMR1", "SB100212TAWMD13VV4TMR2", "SB100212TAWMD14VV4TMR1"))]

# note which OTUs that have less than min.frac
too.small.otus <- rowSums(otu) < min.otu.size

# remove OTUs with small counts
otu <- otu[!too.small.otus, ]

# put the OTU IDs back
sample.ids <- colnames(otu)
otu$OTU_ID <- rownames(otu)
otu <- otu[, c('OTU_ID', sample.ids)]

write.table(otu, "otu.txt", sep="\t", row.names=F, quote=F)
