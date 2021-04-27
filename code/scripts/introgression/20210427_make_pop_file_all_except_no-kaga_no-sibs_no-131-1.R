#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file = args[1]
out_file = args[2]

# Full run

emf_samples = c("hdrr", "hni", "hsok", "javanicus", "melastigma", "ancestor")
vcf_samples = readr::read_tsv(in_file, col_names = "samples")

vcf_samples$population = ifelse(vcf_samples$samples %in% emf_samples,
                            vcf_samples$samples,
                            ifelse(grepl("iCab", vcf_samples$samples),
                                   "icab",
                                   ifelse(grepl("Ho5", vcf_samples$samples),
                                          "ho5",
                                          ifelse(grepl("KW", vcf_samples$samples),
                                                 "kiyosu_wild",
                                                 "mikk"))))

readr::write_tsv(vcf_samples, out_file)
