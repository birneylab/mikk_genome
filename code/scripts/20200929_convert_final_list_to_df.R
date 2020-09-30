#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)

# Collect variables from bash script

file_in <- args[1]
file_out <- args[2]

# Load data
final_list <- readRDS(file_in)

# Reduce down to DF
final_df <- lapply(final_list, function(p1){
  out <- lapply(p1, function(p2){
    # Get D stats
    g_wide_d <- p2[["Genome-wide D"]][["D statistic"]]
    per_chr_d <- sapply(p2[["Per-chromosome"]],
                        function(x) x[["D statistic"]])
    # Get Z scores
    g_wide_z <- p2[["Genome-wide D"]][["Z score"]]
    per_chr_z <- sapply(p2[["Per-chromosome"]],
                        function(x) x[["Z score"]])
    # Get admixture
    admix_f <- p2[["Admixture"]][["f statistic"]]
    per_chr_f <- sapply(p2[["Per-chromosome"]],
                        function(x) x[["Admixture"]][["f statistic"]])
    # Create data frame
    df_out <- data.frame("chr" = c("all", names(per_chr_d)),
                         "d_stat" = c(g_wide_d, per_chr_d),
                         "z_score" = c(g_wide_z, per_chr_z),
                         "admix_f" = c(admix_f, per_chr_f))
    return(df_out)
  })
  # bind rows into single DF
  out <- dplyr::bind_rows(out, .id = "p2")
  return(out)
})

final_df <- dplyr::bind_rows(final_df, .id = "p1")

# tidy up factors, etc for data
final_df$chr <- factor(final_df$chr, levels = c(seq(1, 24), "all"))
## capitalise species/line names
final_df <- final_df %>%
  dplyr::mutate(across(c("p1", "p2"),
                       ~dplyr::if_else(.x %in% c("hni", "hsok"),
                                       str_to_upper(.x),
                                       dplyr::if_else(.x == "hdrr",
                                                      "HdrR",
                                                      .x)))) %>%
  dplyr::rename(P1 = p1,
                P2 = p2)

write.table(final_df, file_out, quote = F, sep = "\t", row.names = F)
