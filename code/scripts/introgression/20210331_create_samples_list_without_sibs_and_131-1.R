#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file = args[1]
out_file = args[2]

# Read in first line of file to get order
samples_orig = scan(in_file, sep = "\t", what = "character", nlines = 1)

# Read in first line of file and split into list
samples = scan(in_file, sep = "\t", what = "character", nlines = 1) %>%
    tibble(SAMPLE = .) %>%
    dplyr::filter(!SAMPLE %in% c("#CHROM", "POS")) %>%
    # Remove iCab, KW and Ho5
    dplyr::filter(grepl("iCab|KW|Ho5", SAMPLE) == F) %>%
    # Remove technical replicates
    dplyr::filter(grepl("-", SAMPLE) == F) %>%
    # Remove 131_1
    dplyr::filter(SAMPLE != "131_1") %>%
    # Split up SAMPLE by LINE and SIB
    dplyr::mutate(LINE = SAMPLE %>% stringr::str_split("_", simplify = T) %>% subset(select = 1),
                  SIB = SAMPLE %>% stringr::str_split("_", simplify = T) %>% subset(select = 2)) %>%
    # Factorise to retain order
    # Split into list
    split(f = .$LINE)

# Pull out sib-lines with the lowest SIB value (e.g. if there are sib lines _2 and _4, take _2)
lapply(samples, function(x){
    if (nrow(x) > 1){
        return(x[which.min(x$SIB), ])
    } else {
        return(x)
    }
}) %>%
    # Bind and write to table
    dplyr::bind_rows() %>%
    dplyr::pull(SAMPLE) %>%
    # Order by original order of samples
    .[order(match(., samples_orig))] %>%
    readr::write_lines(out_file)
