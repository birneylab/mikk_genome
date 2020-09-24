#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

P1 <- args[1]
P2 <- args[2]
dat_file <- args[3]
dir_out <- args[4]

# Source functions

source("mikk_genome/code/scripts/20200916_abbababa_functions.R")

# Read in data

freq_table <- read.table(dat_file,
                         header = T,
                         as.is = T)

# Run abba baba

list_out <- run_abbababa(data = freq_table,
                         P1 = P1,
                         P2 = P2)

# Save list

bname_out <- paste("p1-", P1, "_p2-", P2, ".rds", sep = "")
file_out <- file.path(dir_out, bname_out)
save(list_out, file = file_out)
