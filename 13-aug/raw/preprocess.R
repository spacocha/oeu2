#!/usr/bin/env Rscript

library(dplyr)

min.otu.size <- 0.0025
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
otu <- sum.columns(otu, "M02", c("SB081213TAWMD02VV4TMR1", "SB081213TAWMD02VV4TMR2"))
otu <- rename.columns(otu, "M03", c("SB081213TAWMD03VV4TMR1"))
otu <- rename.columns(otu, "M04", c("SB081213TAWMD04VV4TMR1"))
otu <- sum.columns(otu, "M06", c("SB081213TAWMD06VV4TMR1", "SB081213TAWMD06VV4TMR2"))
otu <- sum.columns(otu, "M07", c("SB081213TAWMD07VV4TMR1", "SB081213TAWMD07VV4TMR2"))
otu <- rename.columns(otu, "M08", c("SB081213TAWMD08VV4TMR1"))
otu <- sum.columns(otu, "M09", c("SB081213TAWMD09VV4TMR1", "SB081213TAWMD09VV4TMR2"))
otu <- sum.columns(otu, "M10", c("SB081213TAWMD10VV4TMR1", "SB081213TAWMD10VV4TMR2"))
otu <- sum.columns(otu, "M11", c("SB081213TAWMD11VV4TMR1", "SB081213TAWMD11VV4TMR2"))
otu <- rename.columns(otu, "M12", c("SB081213TAWMD12VV4TMR1"))
otu <- rename.columns(otu, "M13", c("SB081213TAWMD13VV4TMR1"))
otu <- sum.columns(otu, "M14", c("SB081213TAWMD14VV4TMR1", "SB081213TAWMD14VV4TMR2"))
otu <- sum.columns(otu, "M15", c("SB081213TAWMD15VV4TMR1", "SB081213TAWMD15VV4TMR2"))
otu <- sum.columns(otu, "M17", c("SB081213TAWMD17VV4TMR1", "SB081213TAWMD17VV4TMR2"))
otu <- rename.columns(otu, "M20", c("SB081213TAWMD20VV4TMR1"))
otu <- rename.columns(otu, "M21", c("SB081213TAWMD21VV4TMR1"))
otu <- rename.columns(otu, "M22", c("SB081213TAWMD22VV4TMR1"))
otu <- sum.columns(otu, "MEB", c("SB081213TAWMDEBVV4TMR1", "SB081213TAWMDEBVV4TMR2"))

# normalize by column (i.e., convert to relative abundances)
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

# remove OTUs w/ more than 10% of their reads in the blanks
otu <- otu[otu$MEB  <= max.blank.frac, ]

# remove OTUs w/ more than 5% of their reads in the sediment-y sample (M22)
# otu <- otu[otu$M22 <= max.sediment.frac, ]

# remove extra columns
otu <- otu[, !(colnames(otu) %in% c("MEB"))]

# note which OTUs that have less than min.frac
too.small.otus <- rowSums(otu) < min.otu.size

# remove OTUs with small counts
otu <- otu[!too.small.otus, ]

# put the OTU IDs back
sample.ids <- colnames(otu)
otu$OTU_ID <- rownames(otu)
otu <- otu[, c('OTU_ID', sample.ids)]

write.table(otu, "otu.txt", sep="\t", row.names=F, quote=F)
