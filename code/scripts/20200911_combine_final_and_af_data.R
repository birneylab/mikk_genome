#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

chr <- args[1]
dat_file <- args[2]
af_file <- args[3]
dir_out <- args[4]

# read in data

final_dat <- read.delim(dat_file, header = T)
af_dat <- readr::read_tsv(af_file,
                          col_names = c("chr_pos", "chr", "pos", "ref", "alt", "mikk"),
                          col_types = "ciiccd")

# rename columns

colnames(final_dat)[colnames(final_dat) == "oryzias_latipes"] <- "hdrr" # HdrR

colnames(final_dat) <- gsub("oryzias_latipes_", "", colnames(final_dat)) # HNI and HSOK

colnames(final_dat) <- gsub("oryzias_", "", colnames(final_dat)) # other medaka

# join DFs
joined_dat <- dplyr::inner_join(af_dat,
                                dplyr::select(final_dat,
                                              !c(chr, coord, chr_file)),
                                by = "chr_pos") %>%
              tidyr::drop_na(ancestor) # remove rows without ancestral allele call

# set ancestral and derived alleles and get "frequencies" for each popn

final_dat <- joined_dat %>%
  dplyr::mutate(ancestral = ancestor,
                derived = dplyr::if_else(ancestor == ref,
                                        alt,
                                        ref),
                mikk = dplyr::if_else(ancestral == ref,
                                          mikk,
                                          1 - mikk ),
                mutate(across(hdrr:melastigma,
                              ~dplyr::if_else(.x == derived,
                                              1,
                                              0)))) %>%
  dplyr::select(chr, pos, ancestral, derived,
                mikk, hdrr, hni, hsok, javanicus, melastigma)

# write table

bname_out <- paste(chr, ".txt", sep = "")
file_out <- file.path(dir_out, bname_out)
write.table(final_dat, file_out, quote = F, sep = "\t", row.names = F)
