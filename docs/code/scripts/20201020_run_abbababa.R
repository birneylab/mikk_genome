#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

P1 = args[1]
P2 = args[2]
P3 = args[3]
dat_file = args[4]
dir_out = args[5]

# Source functions

source("mikk_genome/code/scripts/20201020_abbababa_functions.R")

# Read in data

freq_table = read.table(dat_file,
                        header = T,
                        as.is = T,
                        check.names = F)

# Run abba baba

list_out = run_abbababa(data = freq_table,
                        P1 = P1,
                        P2 = P2,
                        P3 = P3)

# Save list

bname_out = paste(P1, "_", P2, "_", P3, ".RData", sep = "")
file_out = file.path(dir_out, bname_out)
saveRDS(list_out, file = file_out)
