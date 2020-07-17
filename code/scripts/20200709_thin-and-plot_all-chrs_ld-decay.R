#!/usr/bin/env Rscript

# Bash script:
# Rscript --vanilla mikk_genome/code/scripts/20200602_ld_decay_plot.R <vcftools.geno.ld file> <output directory> <output prefix (using date today)>

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1]
dir_out <- args[2]
date_today <- args[3]

# Get chromosome number from filename

chr_no <- gsub(".txt", "", basename(file_in))

# Read in data

data <- read.delim(file_in,
                   sep = "\t",
                   header = F,
                   skip = 1)

# Set column names

colnames(data) <- c("chr", "snp_1", "snp_2", "count", "r2")

# Create distance variable

data$distance_kb <- abs(data$snp_1 - data$snp_2) / 1000

# Take subsample of data

set.seed(27)
data <- dplyr::sample_n(data, 5000000)

# Plot

data %>%
  ggplot(aes(distance_kb, r2)) +
  geom_point(size = 0.1, alpha = 0.005) +
  geom_smooth(size = 0.5) +
  theme_bw() +
  facet_wrap(~chr, nrow = 6, ncol = 4) +
  xlab("Distance between SNPs (kb)") +
  ylab(parse(text = "r^2"))

# Save

ggsave(filename = paste(date_today, "_ld_", chr_no, ".png", sep = ""),
       device = "png",
       path = dir_out,
       width = 15,
       height = 15,
       units = "cm",
       dpi = 500)
