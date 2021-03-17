#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)

# Create a parser

p = argparser::arg_parser("Combine ABBA BABA sliding windows data into single file.")

# Add arguments

p = argparser::add_argument(p,
                            arg = "--input",
                            help = "Input files",
                            type = "character",
                            nargs = Inf)

p = argparser::add_argument(p,
                            arg = "--output",
                            help = "Output file",
                            type = "character",
                            nargs = 1)

# Parse args

args = argparser::parse_args(p)

# Get list of files

#files = list.files(args$input,
#                   full.names = T)

files = args$input

# read in files

lapply(files, function(x){
  # read in data
  out_df = read.csv(x)
  # get populations
  split_name = basename(x) %>% 
    stringr::str_split("_", simplify = T)
  p1 = split_name[, 1]
  p2 = split_name[, 2]
  # add to data
  out_df$p1 = p1
  out_df$p2 = p2
  return(out_df)
}) %>% 
  # bind into single DF
  dplyr::bind_rows() %>%
  # write to file
  readr::write_csv(args$output)
