#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

group = args[1]
in_dir = args[2]
out_file = args[3]

# List files

files = list.files(in_dir,
                    pattern = "[0-9]",
                    full.names = T)

# Read in data

data = lapply(files, function(x){
  # read in data
  df = readr::read_tsv(x,
                       col_types = "iiccdiiiii")
  
  return(df)
}) %>% 
    dplyr::bind_rows() %>% # bind into single DF
    dplyr::arrange(chr, pos) # sort by chromosome, then position

# Rename mikk column

mikk_group_name = paste("mikk_", group, sep = "") 
colnames(data)[colnames(data) == "mikk"] = mikk_group_name

# Write file

readr::write_tsv(data, out_file)

