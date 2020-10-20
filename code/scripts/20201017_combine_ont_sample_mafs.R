#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

chr = args[1]
in_dir = args[2]
out_dir = args[3]

# Read files into list

samples = list.files(in_dir) # create vector of samples

df_list = lapply(samples, function(x){ # read into list
  # get path
  target_path = file.path(in_dir, x, paste(chr, ".txt", sep = ""))
  # create names for MAF columns
  suffs = c("", "_a", "_b")
  smp_nms = paste(x, suffs, sep = "")
  # read in table
  smp_df = readr::read_tsv(target_path,
                           col_names = c("chr_pos", "chr", "pos", "ref", "alt", smp_nms),
                           col_types = "ciiccddd")
  return(smp_df)
})
names(df_list) = samples

# Iteratively join data frames by all columns except for MAF columns

df_out = df_list %>%
  purrr::reduce(full_join, by = c("chr_pos", "chr", "pos", "ref", "alt"))

# Write table to file

bname_out = paste(chr, ".txt", sep = "")
file_out = file.path(out_dir, bname_out)
write.table(df_out, file_out, quote = F, sep = "\t", row.names = F)  
