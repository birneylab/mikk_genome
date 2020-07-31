#!/usr/bin/env Rscript

# Bash script:
# Rscript --vanilla mikk_genome/code/scripts/20200715_r2_decay_mean.R <vcftools.geno.ld file> <output directory> <output prefix (using date today)>

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1]
dir_out <- args[2]

# Get chromosome number from filename

chr_no <- gsub(".ld", "", basename(file_in))

# Read in data

data <- readr::read_delim(file_in,
                          delim = " ",
                          col_names = T,
                          trim_ws = T) %>%
            dplyr::select(CHR_A, BP_A, BP_B, R2)

# Set column names

colnames(data) <- c("chr", "snp_1", "snp_2", "r2")

# Create distance variables

data$distance <- abs(data$snp_1 - data$snp_2)
# data$distance_kb <- abs(data$snp_1 - data$snp_2) / 1000


# Create bins
data$bin_dist <- ggplot2::cut_interval(data$distance, n = 500)

# Get left boundary
data$bin <- ggplot2::cut_interval(data$distance, n = 500)
data$bin_bdr <- as.numeric(stringr::str_split(data$bin, ",", simplify = T)[,1 ] %>%
  stringr::str_replace_all("\\(", "") %>%
  stringr::str_replace_all("\\[", ""))

# Group by bin and get mean R^2
out <- data %>%
  dplyr::group_by(bin_bdr) %>%
  dplyr::summarise(mean = mean(r2, na.rm = T))

# Write table
file_out <- file.path(dir_out, paste(chr_no, ".txt", sep = ""))
write.table(x = out,
            file = file_out,
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")
