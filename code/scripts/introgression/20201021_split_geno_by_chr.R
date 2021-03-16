#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

chr = args[1]
in_file = args[2]
out_file = args[3]

# Read in data

df = read.delim(in_file,
                header = F,
                sep = "\t",
                as.is = T)

# Filter for rows matching target chromosome

df = df[df[2] == chr, ]

# Write table

write.table(df, out_file, quote = F, sep = "\t", row.names = F, col.names = F)
