#!/usr/bin/env Rscript

library(dplyr)

# global parameters
otu.count.min <- 2000
otu.sample.min <- 2
k <- 50
prune.min.cor <- 0.75
filter.min.size <- 2
out.fn <- "groups_old.dat"

printf <- function(...) sprintf(...) %>% cat

# read in OTU table
otu <- read.table("otu.txt", header=T, sep="\t", row.names=1)

# remove extra (blanks) columns
otu <- otu[, !(colnames(otu) %in% c("MEB", "MSB"))]

# remove OTUs that have a small number of counts
otu <- otu[rowSums(otu) >= otu.count.min, ]

# remove OTUs that appear in a small number of samples
# (the parens around sum are correct: counting # that are > 0)
# swo> I think this does nothing, all OTUs are present in more than 2
n.samples <- apply(otu, 1, function(x) sum(x > 0))
otu <- otu[n.samples > otu.sample.min, ]

otu.ids <- rownames(otu)

# pool replicate samples, the drop the originals
# swo> this should have been done before counting samples?
otu <- otu %>% rename(M8.1=M8) %>% mutate(M3=M3.1+M3.2, M8=M8.1+M8.2)
otu <- select(otu, -M3.1, -M3.2, -M8.1, -M8.2)

# normalize by column (i.e., convert to relative abundances)
# swo> this should have been done before dropping any rows
otu <- apply(otu, 2, function(x) x / sum(x)) %>% as.data.frame

printf("after preprocessing, there are %d OTUs left\n", dim(otu)[1])

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
      sprintf("group %d has minimum mean correlation %f (i.e., %f), keeping it\n", i, min.mean.cor, 1.0-min.mean.cor) %>% cat
      pruning <- FALSE
    } else {
      # prune the worst performer
      # swo> what about ties? should keep the most abundant
      min.idx <- which(mean.cors == min(mean.cors)) %>% tail(1)
      printf("from group %d, pruning member %d, whose mean was %f (i.e., %f)\n", i, min.idx, min.mean.cor, 1.0-min.mean.cor)
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
out <- data.frame(otu=otu.ids, group=groups)
write.table(out, file=out.fn, quote=FALSE, sep="\t", row.names=F)
