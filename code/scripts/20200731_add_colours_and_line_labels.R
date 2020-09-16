#!/usr/bin/env Rscript

# Bash script:

# Import args from bash script

args <- commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1]

# Read in data

## MIKK
data <- readr::read_delim(file_in,
                              delim = "\t'",
                              trim_ws = T,
                              col_names = F)
