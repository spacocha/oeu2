#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

# global parameters
plot.dir <- "plots"

dat <- read.table("plot.dat", header=T, sep="\t")

for (g in 1:max(dat$group, na.rm=T)) {
  dat %>% filter(group == g) %>%
    ggplot(aes(x=depth, y=value)) + geom_line(aes(group=id)) +
    theme_minimal()
  fn <- sprintf("%s/plot_%02d.pdf", plot.dir, g)
  ggsave(fn)
}
