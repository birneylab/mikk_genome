#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

mikk_all = args[1]
mikk_1 = args[2]
mikk_2 = args[3]
out_file = args[4]

# Read in data

mikk_all = readr::read_tsv(mikk_all,
                           col_types = "iiccdiiiii")

mikk_1 = readr::read_tsv(mikk_1,
                           col_types = "iiccdiiiii") %>%
                dplyr::select(chr, pos, mikk_1)

mikk_2 = readr::read_tsv(mikk_2,
                           col_types = "iiccdiiiii") %>%
                dplyr::select(chr, pos, mikk_2)

# Bind all

final_df = mikk_all %>%
    dplyr::full_join(mikk_1, by = c("chr", "pos")) %>%
    dplyr::full_join(mikk_2, by = c("chr", "pos")) %>%
    dplyr::select(chr, pos, ancestral, derived, mikk = mikk_all, mikk_a = mikk_1, mikk_b = mikk_2, everything())

# Write to file

readr::write_tsv(final_df, out_file)