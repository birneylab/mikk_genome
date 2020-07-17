#!/usr/bin/env Rscript

# Bash script:
# Rscript --vanilla mikk_genome/code/scripts/20200715_r2_decay_mean.R <vcftools.geno.ld file> <output directory> <output prefix (using date today)>

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1] # testing with "ld/20200707_panel_maf-0.10_window-50kb_no-missing_maf-max-0.90/18.geno.ld", which has the most lines at 82M
dir_out <- args[2] # testing with

# Get chromosome number from filename

chr_no <- gsub(".geno.ld", "", basename(file_in))

# Read in data

data <- read.delim(file_in,
                   sep = "\t",
                   header = F,
                   skip = 1)

# Set column names

colnames(data) <- c("chr", "snp_1", "snp_2", "count", "r2")

# Create distance variables

data$distance <- abs(data$snp_1 - data$snp_2)
# data$distance_kb <- abs(data$snp_1 - data$snp_2) / 1000

# Filter for SNPs <= 1000 bp away from each other
data <- data %>%
  dplyr::filter(distance <= 1000)

# Create bins
data$bin_dist <- ggplot2::cut_interval(data$distance, n = 50)

# Get left boundary
data$bin <- ggplot2::cut_interval(data$distance, n = 50)
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
