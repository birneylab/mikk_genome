#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file = args[1]
target_mikk = args[2]
out_file = args[3]

# Read in file, exclude target_mikk line, then write to file

data = readr::read_lines(in_file)

data = data[which(data != target_mikk)]

readr::write_lines(data, out_file)

