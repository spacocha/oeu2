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
otu <- rename.columns(otu, "M00", c("X813080Q1"))
otu <- rename.columns(otu, "M01", c("X8130815Q1"))
otu <- sum.columns(otu, "M03", c("X813083Q7", "X813083Q8"))
otu <- rename.columns(otu, "M04", c("X813084Q1"))
otu <- rename.columns(otu, "M05", c("X813085Q1"))
otu <- rename.columns(otu, "M06", c("X813086Q1"))
otu <- rename.columns(otu, "M07", c("X813087Q1"))
otu <- sum.columns(otu, "M08", c("X813088Q1", "X813088Q2"))
otu <- rename.columns(otu, "M09", c("X813089Q2"))
otu <- rename.columns(otu, "M10", c("X8130810Q1"))
otu <- rename.columns(otu, "M11", c("X8130811Q1"))
otu <- rename.columns(otu, "M12", c("X8130812Q1"))
otu <- rename.columns(otu, "M13", c("X8130813Q1"))
otu <- rename.columns(otu, "M14", c("X8130814Q1"))
otu <- rename.columns(otu, "M15", c("X8130815Q1_1"))
otu <- rename.columns(otu, "M16", c("X8130816Q1"))
otu <- rename.columns(otu, "M17", c("X8130817Q3"))
otu <- rename.columns(otu, "M19", c("X8130819Q3"))
otu <- rename.columns(otu, "M20", c("X8130820Q1"))
otu <- rename.columns(otu, "M21", c("X8130821Q1"))
otu <- rename.columns(otu, "M22", c("X8130822Q1"))
otu <- sum.columns(otu, "MEB", c("X81308SBQ1", "X81308EBQ1"))

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
