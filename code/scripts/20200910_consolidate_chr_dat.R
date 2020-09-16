#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

dir_in <- args[1]
file_out <- args[2]

# Get list of files

files <- list.files(dir_in, full.names = T)

# read data files into list and process

dat_list <- lapply(files, function(x){
  # read in data
  df <- read.delim(x, header = T)
  # change `chr` and `coord` columns into integers
  df$chr <- as.integer(df$chr)
  df$coord <- as.integer(df$coord)
  # change name of ancestor in header
  colnames(df)[grep("Ancestor", colnames(df))] <- "ancestor"
  # capitalise all letters
  df <- df %>%
    dplyr::mutate(across(.cols = c(-chr, -coord),
                         .fns = toupper))
  return(df)
})

# tidy data

dat_df <- dplyr::bind_rows(dat_list) %>% # bind into single DF
  dplyr::arrange(coord) %>%  # sort by coordinate
  tidyr::unite("chr_pos", chr, coord, sep = ":", remove = F) # create new column

# write table

write.table(dat_df, file_out, quote = F, sep = "\t", row.names = F)
