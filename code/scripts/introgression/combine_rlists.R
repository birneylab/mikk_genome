#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

in_dir = args[1]
out_file = args[2]

# get list of target files
files <- list.files(in_dir,
                    full.names = T)

# read into one list
final_lst = lapply(files, readRDS)

# Reduce down to DF
final_df = lapply(final_lst, function(file){
  # Get D stats
  g_wide_d <- file[["Genome-wide D"]][["D statistic"]]
  per_chr_d <- sapply(file[["Per-chromosome"]],
                      function(x) x[["D statistic"]])
  # Get Z scores
  g_wide_z <- file[["Genome-wide D"]][["Z score"]]
  per_chr_z <- sapply(file[["Per-chromosome"]],
                      function(x) x[["Z score"]])
  # Get admixture
  admix_f <- file[["Admixture"]][["f statistic"]]
  per_chr_f <- sapply(file[["Per-chromosome"]],
                      function(x) x[["Admixture"]][["f statistic"]])
  # Get confidence intervals
  g_wide_ci_lower <- file[["Admixture"]][["Confidence interval"]][["lower"]]
  g_wide_ci_upper <- file[["Admixture"]][["Confidence interval"]][["upper"]]
  per_chr_ci_lower <- sapply(file[["Per-chromosome"]],
                             function(x) x[["Admixture"]][["Confidence interval"]][["lower"]])
  per_chr_ci_upper <- sapply(file[["Per-chromosome"]],
                             function(x) x[["Admixture"]][["Confidence interval"]][["upper"]])

  # Create data frame
  df_out <- data.frame("p1" = file[["Populations"]][["P1"]],
                       "p2" = file[["Populations"]][["P2"]],
                       "p3" = file[["Populations"]][["P3"]],
                       "chr" = c("all", names(per_chr_d)),
                       "d_stat" = c(g_wide_d, per_chr_d),
                       "z_score" = c(g_wide_z, per_chr_z),
                       "admix_f" = c(admix_f, per_chr_f),
                       "f_ci_lower" = c(g_wide_ci_lower, per_chr_ci_lower),
                       "f_ci_upper" = c(g_wide_ci_upper, per_chr_ci_upper),
                       stringsAsFactors = F)
  return(df_out)
}) %>% 
  dplyr::bind_rows()

# capitalise species / line names / colnames
rename_vec = c("HdrR", "HNI", "HO5", "HSOK", "iCab", "KW", "MIKK")
names(rename_vec) = c("hdrr", "hni", "ho5", "hsok", "icab", "kiyosu_wild", "mikk")

final_df = final_df %>% 
  dplyr::mutate(dplyr::across(c("p2", "p3"),
                              ~dplyr::recode(.x, !!!rename_vec)))

colnames(final_df) = dplyr::recode(colnames(final_df),
                                   p1 = "P1",
                                   p2 = "P2",
                                   p3 = "P3")

# save
readr::write_tsv(final_df, out_file)

