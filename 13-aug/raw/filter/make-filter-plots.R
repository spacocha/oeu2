#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)
library(ggplot2)

theme_set(theme_minimal())

min.otu.size <- 100

otu <- read.table("../otu_rename.txt", header=T, sep="\t", row.names=1)

dat <- data.frame(total.counts=apply(otu, 1, sum),
                  sediment.counts=otu$M22) %>%
  filter(total.counts > 0) %>%
  mutate(frac.sediment=sediment.counts/total.counts)

# make the sediment fraction plot
p <- dat %>% ggplot(aes(x=total.counts, y=frac.sediment)) +
  geom_point(size=1) +
  scale_x_log10() +
  scale_y_continuous(breaks=c(0, 0.05, 0.5, 1.0)) +
  xlab("counts in all samples") + ylab("fraction of counts in near-bottom sample") +
  geom_segment(aes(x=min.otu.size, xend=min.otu.size, y=-0.02, yend=0.05), col='red') +
  geom_segment(aes(x=min.otu.size, xend=max(dat$total.counts), y=0.05, yend=0.05), col='red')

ggsave('sediment-filter.pdf', p)
