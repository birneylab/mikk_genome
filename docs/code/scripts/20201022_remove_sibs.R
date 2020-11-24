#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file = args[1]
sibs_file = args[2]
out_file = args[3]

# Read in files

df = read.table(in_file,
                header = T,
                sep = "\t",
                as.is = T,
                comment.char = "*",
                check.names = F)

# Remove siblings

## Read in siblings file
excl_sibs = scan("mikk_genome/data/20200227_panel_lines_excluded.txt", character())
## Substitute dashes for underscores
excl_sibs = gsub("-", "_", excl_sibs)
## Remove siblings
df_nosibs = df %>%
  dplyr::select(-all_of(excl_sibs))

# Write table

write.table(df_nosibs, out_file, quote = F, sep = "\t", row.names = F)
