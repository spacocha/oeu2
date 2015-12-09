#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

# global parameters
min.otu.size <- 1000
prune.min.cor <- 0.75
filter.min.size <- 2
max.blank.frac <- 0.10
max.sediment.frac <- 0.05
top.otus <- 400

# read in OTU table
otu <- read.table("../otu.txt", header=T, sep="\t", row.names=1)

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

# sort OTUs by abundance (by summing across rows)
otu <- otu[order(rowSums(otu), decreasing=TRUE), ]

# swo> there are 484 OTUs now
# take the most abundant OTUs only
# otherwise, you get a funny drop-off in r.a. at ~475
otu <- otu[1:top.otus, ]

depths <- sapply(names(otu), function(x) as.numeric(gsub("M", "", x)))
mean.depths <- apply(otu, 1, function(x) weighted.mean(depths, w=x))

# normalize the OTUs by row
otu <- apply(otu, 1, function(x) x / sum(x)) %>% t %>% as.data.frame

# compute the dissimilarities
euc.dist <- dist(otu, method="euclidean")
cor.dist <- cor(t(otu))

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

cluster.at <- function(k) {
  # do the clustering
  fit <- hclust(euc.dist, method='ward.D')
  groups <- cutree(fit, k=k)

  # prune each group
  for (i in 1:k) groups <- prune(groups, i)

  # filter the groups (i.e., remove ones below a certain abundance threshold)
  for (i in 1:k) {
    group.size <- tabulate(groups)[i]
    if (group.size < filter.min.size) {
      groups[groups == i] <- NA
    }
  }

  # refactor the groups
  while (any(tabulate(groups) == 0)) {
    unused.val <- which(tabulate(groups) == 0) %>% head(1)
    idx <- groups > unused.val & is.finite(groups)
    groups[idx] <- groups[idx] - 1
  }

  # create the output information
  n.groups <- max(groups, na.rm=T)
  mean.group.size <- mean(tabulate(groups))
  mean.var.mean <- sapply(1:max(groups, na.rm=T), function(g) var(mean.depths[is.finite(groups) & groups == g])) %>% mean
  total.otus <- sum(!is.na(groups))

  # return a vector of this info
  c(k, n.groups, mean.group.size, mean.var.mean, total.otus)
}

cak <- sapply(5:100, cluster.at) %>% t %>% as.data.frame
names(cak) <- c("k", "n.oeus", "mean.size", "mvm", "total.otus")
write.table(cak, file="cak.dat", sep="\t", quote=F, row.names=F)

ggplot(cak, aes(x=k, y=n.oeus)) + geom_point() + theme_minimal() + ylab("number of final OEUs")
ggsave("plots/noeus.pdf")
ggplot(cak, aes(x=k, y=mean.size)) + geom_point() + theme_minimal() + ylab("mean OTUs per OEU")
ggsave("plots/meansize.pdf")
ggplot(cak, aes(x=k, y=mvm)) + geom_point() + theme_minimal() + ylab("mean variance in mean depth")
ggsave("plots/mvm.pdf")
ggplot(cak, aes(x=k, y=total.otus)) + geom_point() + theme_minimal() + ylab("total OTUs in all OEUs")
ggsave("plots/totalotus.pdf")
