#!/usr/bin/env Rscript

library(dplyr)

# global parameters
k <- 50
prune.min.cor <- 0.75
filter.min.size <- 2
out.fn <- "groups.dat"

# read in OTU table
otu <- read.table("otu.txt", header=T, sep="\t", row.names=1)

# pool replicate samples, the drop the originals
otu$M3.pool <- otu$M3.1 + otu$M3.2
otu$M8.pool <- otu$M8 + otu$M8.2
otu <- otu[, !(colnames(otu) %in% c("M3.1", "M3.2", "M8", "M8.2"))]

# note which OTUs that have fewer than 1000 counts
too.small.otus <- rowSums(otu) < 1000

# normalize by column (i.e., convert to relative abundances)
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

# remove OTUs with small counts
otu <- otu[!too.small.otus, ]

# remove OTUs w/ more than 10% of their reads in the blanks
otu <- otu[otu$MEB + otu$MSB <= 0.10, ]

# remove OTUs w/ more than 5% of their reads in the sediment-y sample (M22)
otu <- otu[otu$M22 <= 0.05, ]

# remove extra columns
otu <- otu[, !(colnames(otu) %in% c("MEB", "MSB", "M22"))]

# sort OTUs by abundance (by summing across rows)
otu <- otu[order(rowSums(otu), decreasing=TRUE), ]

# swo> there are 484 OTUs now
# take the most abundant 300 OTUs only
# otherwise, you get a funny drop-off in r.a. at ~475
otu <- otu[1:300, ]

# normalize the OTUs by row
otu <- apply(otu, 1, function(x) x / sum(x)) %>% t %>% as.data.frame

# compute the dissimilarities
euc.dist <- dist(otu, method="euclidean")
cor.dist <- cor(t(otu))

# do the clustering
fit <- hclust(euc.dist, method='ward.D')
groups <- cutree(fit, k=k)

# how-to prune each group
prune <- function(groups, i) {
  # prune the i-th group
  pruning <- TRUE
  while (pruning) {
    members <- which(groups == i)
    if (length(members) == 1) break
    idx <- as.numeric(members)
    
    # for each member of the group, find its mean correlation
    these.cors <- cor.dist[idx, idx]
    diag(these.cors) <- NA
    mean.cors <- rowMeans(these.cors, na.rm=T)
    
    # if the minimum mean correlation is above a threshold, we're done
    min.mean.cor <- min(mean.cors)
    if (min.mean.cor > prune.min.cor) {
      sprintf("group %d has minimum mean correlation %f, keeping it\n", i, min.mean.cor) %>% cat
      pruning <- FALSE
    } else {
      # prune the worst performer
      # swo> what about ties? should keep the most abundant
      min.idx <- which(mean.cors == min(mean.cors)) %>% tail(1)
      groups[members[min.idx]] <- NA
    }
  }
  # return the new groups vector, which may have had some elements replaced by NAs
  groups
}

# prune each group
for (i in 1:k) groups <- prune(groups, i)

# filter the groups (i.e., remove ones below a certain abundance threshold)
for (i in 1:k) {
  group.size <- tabulate(groups)[i]
  if (group.size < filter.min.size) {
    sprintf("filtering group %d, which has size %d\n", i, group.size) %>% cat
    groups[groups == i] <- NA
  }
}

# refactor the groups
while (any(tabulate(groups) == 0)) {
  unused.val <- which(tabulate(groups) == 0) %>% head(1)
  sprintf("refactoring to reuse empty name %d\n", unused.val) %>% cat
  idx <- groups > unused.val & is.finite(groups)
  groups[idx] <- groups[idx] - 1
}

# prepare an output table
out <- data.frame(otu=names(groups), group=groups)
write.table(out, file=out.fn, quote=FALSE, sep="\t", row.names=F)
