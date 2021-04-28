#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

emf_in_file = args[1]
vcf_gen_file = args[2]
vcf_samples = args[3]
out_file = args[4]

# Create colnames vector for vcf

samples = scan(vcf_samples, character())
vcf_cols = c("chr_pos", "chr", "pos", "ref", "alt", samples)

# Read in files

emf_df = read.table(emf_in_file,
                    header = T,
                    sep = "\t",
                    as.is = T)

vcf_df = read.table(vcf_gen_file,
                    header = F,
                    sep = "\t",
                    as.is = T,
                    col.names = vcf_cols,
                    check.names = F) %>%
  dplyr::select(-c(ref, alt))

# Make recode vector

recode_vec = c("A/A", "G/G", "T/T", "C/C")
names(recode_vec) = c("A", "G", "T", "C")

# Select relevant cols, remove NAs and recode to geno format

final_df = emf_df %>%
  # select relevant cols
  dplyr::select(chr_pos,
                chr,
                pos = coord,
                hdrr = oryzias_latipes,
                hni = oryzias_latipes_hni,
                hsok = oryzias_latipes_hsok,
                ancestor) %>%
  # replace NA with "N/N"
  dplyr::mutate(across(hdrr:ancestor,
                ~replace_na(.x, "N/N"))) %>%
  # recode haploid as diploid in .geno format
  dplyr::mutate(across(hdrr:ancestor,
                ~dplyr::recode(.x, !!!recode_vec))) %>%
  # join with VCF file
  dplyr::inner_join(vcf_df,
                    by = c("chr_pos", "chr", "pos")) %>%
  # remove chr_pos column and rename chr and pos columns
  dplyr::select("#CHROM" = chr,
                "POS" = pos,
                everything()) %>%
  dplyr::select(!chr_pos)

# Write table

write.table(final_df, out_file, quote = F, sep = "\t", row.names = F)
