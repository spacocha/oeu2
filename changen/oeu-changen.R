#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

min.k <- 5
max.k <- 100

mean.depths <- read.table("mean-depths.txt", header=T, sep="\t")$mean.depth

at.k <- function(k) {
  res <- system2("../oeu.R", args=c("-k", as.character(k), "-n", "1000", "otu.txt", "groups.dat"))
  groups <- read.table("groups.dat", header=T, sep="\t")$group
  n.groups <- max(groups, na.rm=T)
  n.otus <- sum(!is.na(groups))
  avg.size <- tabulate(groups) %>% mean
  
  # variance in the mean of depths
  vmd <- function(g) var(mean.depths[groups == g], na.rm=T)
  
  # mean variance in the mean of depths
  mvm <- sapply(1:n.groups, vmd) %>% mean
  
  c(k, n.groups, n.otus, avg.size, mvm)
}

out <- sapply(min.k:max.k, at.k) %>% t %>% as.data.frame
colnames(out) <- c("k", "n.oeus", "n.otus", "avg.size", "mvm")

# save the output
write.table(out, "res.dat", sep="\t", row.names=F, quote=F)

# make some simple plots
theme_set(theme_minimal())
p <- ggplot(out, aes(x=k, y=n.oeus)) +
  geom_point() +
  xlab("number of initial candidate OEUs") +
  ylab("number of final OEUs")
ggsave("plot/noeus.pdf", p)

p <- ggplot(out, aes(x=k, y=n.otus)) + geom_point() +
  xlab("number of initial candidate OEUs") +
  ylab("number of final OTUs in all OEUs")
ggsave("plot/notus.pdf", p)

p <- ggplot(out, aes(x=k, y=avg.size)) + geom_point() +
  xlab("number of initial candidate OEUs") +
  ylab("mean OTUs per OEU")
ggsave("plot/avgsize.pdf", p)

p <- ggplot(out, aes(x=k, y=mvm)) + geom_point() +
  xlab("number of intial candidate OEUs") +
  ylab("mean variance of mean depth within OEUs")
ggsave("plot/mvm.pdf", p)
