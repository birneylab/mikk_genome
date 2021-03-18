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

# Read files 

files = args$input

# Process and write file

lapply(files, function(x){
  out = read.table(x,
                   header = T, 
                   as.is = T, 
                   comment.char = "*",
                   check.names = F)
  return(out)
}) %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(chr = "#CHROM",
                  pos = "POS",
                  everything()) %>% 
    dplyr::mutate(dplyr::across(hni:hsok,
                                ~ifelse(.x == "N/N", NA, 1))) %>%
    readr::write_tsv(args$output)
