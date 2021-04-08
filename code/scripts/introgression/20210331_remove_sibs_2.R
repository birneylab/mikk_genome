#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file = args[1]
samples_keep_file = args[2]
out_file = args[3]

# Read in files

df = read.table(in_file,
                header = T,
                sep = "\t",
                as.is = T,
                comment.char = "*",
                check.names = F)

# Read in samples to keep

samples_keep = scan(samples_keep_file, what = character())

# Remove them

df_nosibs = df %>%
  dplyr::select("#CHROM", POS, all_of(samples_keep))

# Write table

write.table(df_nosibs, out_file, quote = F, sep = "\t", row.names = F)
