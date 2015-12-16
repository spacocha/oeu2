#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))

p <- OptionParser(usage="%prog [options] input_otu_table output_groups")
p <- add_option(p, c("-k", "--candidates"), type="integer", default=50, 
  help="number of initial/candidate OEUs (default: %default)", metavar="int")
p <- add_option(p, c("-c", "--correlation"), type="numeric", default=0.75,
  help="minimum average correlation before pruning an OTU from an OEU (default: %default)", metavar="float")
p <- add_option(p, c("-n", "--n_otus"), type="integer", default=400,
  help="number of OTUs to include in the analysis (default: %default)", metavar="int")
p <- add_option(p, c("-f", "--filter"), type="integer", default=2,
  help="minimum number of OTUs in an OEU (default: %default)", metavar="int")
args <- parse_args(p, positional_arguments=2)
opt <- args$options
in.fn <- args$args[1]
out.fn <- args$args[2]

# global parameters
k <- opt$candidates
prune.min.cor <- opt$correlation
filter.min.size <- opt$filter
top.otus <- opt$n_otus

otu <- read.table(in.fn, header=T, sep="\t", row.names=1)

# sort OTUs by abundance (by summing across rows)
otu <- otu[order(rowSums(otu), decreasing=TRUE), ]

# take the most abundant OTUs only
n.otus <- dim(otu)[1]
if (n.otus < top.otus) {
  sprintf("only %d OTUs left after filtering, fewer than max %d, proceeding anyway", n.otus, top.otus) %>% cat
} else {
  otu <- otu[1:top.otus, ]
}

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
      # prune the worst performer, keeping the most abundant if there's a tie
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
