#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1]
file_out <- args[2]

# Read in data

df <- readr::read_tsv(file_in,
                      col_names = c("chr", "pos", "ref", "alt", "mikk"),
                      col_types = "ciccd")

# Create chr_pos column

df$chr_pos <- paste(df$chr, df$pos, sep = ":")

# Create pop_a and pop_b columns by randomly assigning 0 or 1 when heterozygous

df$mikk_a <- ifelse(df$mikk == 0.5,
                    sample(c(0,1), size = nrow(df[df$mikk == 0.5, ]), replace = T),
                    df$mikk)

df$mikk_b <- ifelse(df$mikk != 0.5,
                    df$mikk,
                    ifelse(df$mikk_a == 0,
                           1,
                           0))

# Reorder

df <- df %>%
  dplyr::select(chr_pos, everything())

# Write table
write.table(df, file_out, quote = F, sep = "\t", row.names = F, col.names = F)
