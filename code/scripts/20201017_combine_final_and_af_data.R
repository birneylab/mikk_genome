#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

chr = args[1]
in_file_mik = args[2]
in_file_sub = args[3]
in_file_ont = args[4]
in_file_emf = args[5]
out_dir = args[6]

# Read in and combine MIKK, subpop and ONT files

## Create vector of files
target_maf_files = c(in_file_mik, in_file_sub, in_file_ont)
## read in 3 files
maf_df = lapply(target_maf_files, function(x){
  out_df = read.table(x, header = T, as.is = T, check.names = F)
  return(out_df)
})
## join into single df
maf_df = maf_df %>%
  purrr::reduce(full_join, by = c("chr_pos", "chr", "pos", "ref", "alt"))

# Read in EMF data

emf_df = read.delim(in_file_emf, header = T, as.is = T)
## rename columns
colnames(emf_df)[colnames(emf_df) == "oryzias_latipes"] = "hdrr" # HdrR
colnames(emf_df) = gsub("oryzias_latipes_", "", colnames(emf_df)) # HNI and HSOK
colnames(emf_df) = gsub("oryzias_", "", colnames(emf_df)) # other medaka

# Join DFs
joined_dat = dplyr::inner_join(dplyr::select(emf_df,
                                             !c(chr, coord, chr_file)),
                               maf_df,
                               by = "chr_pos")

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
joined_dat = joined_dat %>%
  dplyr::mutate(across(hdrr:melastigma,
                ~ifelse(.x == ref | .x == alt,
                .x,
                NA))) %>%
  dplyr::mutate(across(mikk:last_col(), # make all columns from 'mikk' numeric
                       as.numeric)) %>%
  tidyr::drop_na(ancestor) # and remove rows with NA in the 'ancestor' column

# set ancestral and derived alleles and get "frequencies" for each popn
## get name of final column
final_col = tail(colnames(joined_dat), 1)
## run over df
final_dat = joined_dat %>%
  dplyr::mutate(ancestral = ancestor,
                derived = dplyr::if_else(ancestral == ref,
                                         alt,
                                         ref),
                mutate(across(hdrr:melastigma,
                             ~dplyr::if_else(.x == derived,
                                             1,
                                             0))),
                mutate(across(mikk:all_of(final_col),
                              ~dplyr::if_else(ancestral == ref,
                                              .x,
                                              1 - .x )))) %>%
  dplyr::select(chr, pos, ancestral, derived,
                hdrr, hni, hsok, javanicus, melastigma,
                everything(), -c(chr_pos, ancestor, ref, alt))

# write table

bname_out = paste(chr, ".txt", sep = "")
file_out = file.path(out_dir, bname_out)
write.table(final_dat, file_out, quote = F, sep = "\t", row.names = F)
