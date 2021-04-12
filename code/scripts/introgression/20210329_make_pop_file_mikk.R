#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(argparser)
library(tidyverse)

# Collect variables from bash script

# Create a parser
p = argparser::arg_parser("Create .pop.txt file for ABBABABAsliding.py")

# Add arguments
p = argparser::add_argument(p,
                            arg = "--input",
                            help = "Input .pop.txt file",
                            type = "character",
                            nargs = 1)

p = argparser::add_argument(p,
                            arg = "--abba",
                            help = "Lines always to be included for ABBA-BABA",
                            type = "character",
                            nargs = Inf)

p = argparser::add_argument(p,
                            arg = "--mikk",
                            help = "MIKK lines to include for ABBA-BABA",
                            type = "character",
                            nargs = Inf)

p = argparser::add_argument(p,
                            arg = "--output",
                            help = "Output file",
                            type = "character",
                            nargs = 1)

# Parse args
args = argparser::parse_args(p)

# Get target columns
target_cols = c(args$abba, args$mikk)

# Read in file
readr::read_tsv(args$input) %>%
       dplyr::filter(samples %in% target_cols) %>%
       # Substitute dash for underscore
       dplyr::mutate(samples = stringr::str_replace(samples, "-", "_")) %>%
       readr::write_tsv(args$output)