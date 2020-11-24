#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

chr <- args[1]
dat_file <- args[2]
af_file_mikk <- args[3]
af_file_ont <- args[4]
af_file_subpop <- args[5]
dir_out <- args[6]

# read in data

final_dat <- read.delim(dat_file, header = T, as.is = T)
af_dat <- readr::read_tsv(af_file,
                          col_names = c("chr_pos", "chr", "pos", "ref", "alt", "mikk", "mikk_a", "mikk_b"),
                          col_types = "ciiccddd")

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

# replace non-ref or alt alleles
### check how many alleles in MAF files
#allele_count = apply(joined_dat[, c("ref", "alt", "hdrr", "hni", "hsok", "ancestor", "javanicus", "melastigma")], 1, function(x){
#  length(unique(na.omit(x)))
#})
### get indices of rows with more than 2 alleles
#mult_alleles = which(allele_count > 2)
### how many?
#length(mult_alleles) # 114436 / 482881 = ~23.6% of chr 10 file
## replace non-ref or -alt alleles with NA
joined_dat <- joined_dat %>%
  dplyr::mutate(across(hdrr:melastigma,
                ~ifelse(.x == ref | .x == alt,
                .x,
                NA))) %>%
  tidyr::drop_na(ancestor) # and remove rows with NA in the 'ancestor' column

# set ancestral and derived alleles and get "frequencies" for each popn

final_dat <- joined_dat %>%
  dplyr::mutate(ancestral = ancestor,
                derived = dplyr::if_else(ancestral == ref,
                                         alt,
                                         ref),
                mutate(across(mikk:mikk_b,
                              ~dplyr::if_else(ancestral == ref,
                                      .x,
                                      1 - .x ))),
                mutate(across(hdrr:melastigma,
                              ~dplyr::if_else(.x == derived,
                                              1,
                                              0)))) %>%
  dplyr::select(chr, pos, strand, ancestral, derived,
                mikk, mikk_a, mikk_b,
                hdrr, hni, hsok, javanicus, melastigma)

# write table

bname_out <- paste(chr, ".txt", sep = "")
file_out <- file.path(dir_out, bname_out)
write.table(final_dat, file_out, quote = F, sep = "\t", row.names = F)
