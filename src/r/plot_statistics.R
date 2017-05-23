#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# read data
raw_data <- fread("output/assembly_statistics/stats.txt")

ggplot(raw_data, aes(x = basename(filename), y = scaf_bp)) +
    theme(axis.text.x = element_text(angle = 30)) +
    geom_col()
