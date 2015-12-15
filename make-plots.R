#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(reshape2)

# global parameters
groups.fn <- "groups.dat"
otu.fn <- "08_aug/otu.txt"
plot.dir <- "plots"

# read in the groups data and OTU table
group <- read.table(groups.fn, header=T, sep="\t")
otu <- read.table(otu.fn, header=T, sep="\t", row.names=1)

# rescale the rows (i.e., make profiles)
otu <- apply(otu, 1, function(x) x / sum(x)) %>% t %>% as.data.frame
otu$otu <- rownames(otu)

# create a data frame with entries otu, group, sample, value (= ra)
dat <- merge(group, otu) %>% melt(id.vars=c("otu", "group"), variable.name="sample")

# add the depths by parsing the sample names
dat$depth <- sapply(dat$sample %>% as.list, function(x) sub("M", "", x) %>% as.numeric)

write.table(dat, 'tmp2', row.names=F, quote=F, sep="\t")

# plot each piece of the output
for (g in 1:max(dat$group, na.rm=T)) {
  foo <- filter(dat, group == g)
  sprintf("group %d has %d rows\n", g, dim(foo)[1]) %>% cat
  dat %>% filter(group == g) %>%
    ggplot(aes(x=depth, y=value)) + geom_line(aes(group=otu)) +
    theme_minimal()
  fn <- sprintf("%s/plot_%02d.pdf", plot.dir, g)
  ggsave(fn)
}
