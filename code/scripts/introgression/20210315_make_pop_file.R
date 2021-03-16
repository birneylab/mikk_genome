#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_file_full_run = args[1]
in_file_no_sibs = args[2]
out_file_full_run = args[3]
out_file_no_sibs = args[4]

# Full run

emf_samples = c("hdrr", "hni", "hsok", "javanicus", "melastigma", "ancestor")
vcf_samples = scan(in_file_full_run, character())

pop_out = data.frame("samples" = c(emf_samples, vcf_samples),stringsAsFactors = F)

pop_out$population = ifelse(pop_out$samples %in% emf_samples,
                            pop_out$samples,
                            ifelse(grepl("iCab", pop_out$samples),
                                   "icab",
                                   ifelse(grepl("Ho5", pop_out$samples),
                                          "ho5",
                                          ifelse(grepl("KW", pop_out$samples),
                                                 "kiyosu_wild",
                                                 "mikk"))))

readr::write_tsv(pop_out, out_file_full_run)

# No sibs

no_sibs = scan(in_file_no_sibs, character())
pop_out_no_sibs = pop_out %>% # filter
       dplyr::filter(samples %in% no_sibs | samples %in% emf_samples)

readr::write_tsv(pop_out_no_sibs, out_file_no_sibs)