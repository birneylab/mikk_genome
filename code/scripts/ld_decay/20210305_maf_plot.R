#!/usr/bin/env Rscript

# Bash script:
# Rscript --vanilla mikk_genome/code/scripts/20200727_maf_plot.R maf/20200727_mikk_no-missing.frq maf/20200727_1kg_no-missing.frq plots 20200727_maf_MIKKv1KG

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)
library(cowplot)

# Collect variables from bash script

in_mikk <- args[1]
in_1kg <- args[2]
out_file <- args[3]

# Read in data

## MIKK
maf_mikk <- readr::read_delim(in_mikk,
                             delim = " ",
                             trim_ws = T,
                             col_types = cols_only(MAF = col_double()))
maf_mikk$dataset <- "MIKK"

## 1KG
maf_1kg <- readr::read_delim(in_1kg,
                             delim = " ",
                             trim_ws = T,
                             col_types = cols_only(MAF = col_double()))
maf_1kg$dataset <- "1KG"

## Bind
maf_final <- rbind(maf_mikk, maf_1kg)

# Plot
maf_final %>%
  ggplot() +
    geom_histogram(aes(x = MAF,
                       y=0.01*..density..,
                       fill = dataset),
                   binwidth = 0.01) +
    theme_cowplot() +
    guides(fill = F) +
    facet_wrap(~dataset, nrow = 1, ncol = 2) +
    xlab("Minor allele frequencies") +
    ylab("Density") +
    theme(strip.background = element_blank()) #+
#  scale_fill_manual(values = c(`1KG` = "#FC4E07",
#                               MIKK = "#360568"))

# Save
ggsave(filename = paste(out_file, ".svg", sep = ""),
       device = "svg",
       width = 21.75,
       height = 8,
       units = "cm")
