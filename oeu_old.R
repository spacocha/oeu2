#!/usr/bin/env Rscript

library(dplyr)

# global parameters
otu.count.min <- 2000
otu.sample.min <- 2
k <- 50
prune.max.div <- 0.25
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

# normalize the OTUs by row
otu <- apply(otu, 1, function(x) x / sum(x)) %>% t %>% as.data.frame

# compute the dissimilarities: euclidean distance and correlation "divergence"
# (i.e., 1 - correlation)
euc.dist <- dist(otu, method="euclidean")
cor.div <- 1.0 - cor(t(otu))

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
    
    # for each member of the group, find its mean correlation divergence
    these.divs <- cor.div[idx, idx]
    diag(these.divs) <- NA
    mean.divs <- rowMeans(these.divs, na.rm=T)
    
    # if the maximum mean correlation is below the threshold, we're done
    max.mean.div <- max(mean.divs)
    if (max.mean.div < prune.max.div) {
      printf("group %d has max mean divergence %f, keeping it\n", i, max.mean.div)
      pruning <- FALSE
    } else {
      # prune the worst performer
      # swo> what about ties? should keep the most abundant, but following the old
      # standard of numpy.argmax, which grabs FIRST index
      max.idx <- which(mean.divs == max(mean.divs)) %>% head(1)
      printf("from group %d, pruning member %d (id %d), whose mean was %f\n", i, max.idx, members[max.idx], max.mean.div)
      groups[members[max.idx]] <- NA
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
