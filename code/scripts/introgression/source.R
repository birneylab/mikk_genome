################################
# Libraries
################################

library(renv)
library(ape)
library(geiger)
library(tidyverse)
library(cowplot)
library(karyoploteR)
library(GenomicRanges)
library(biomaRt)
library(circlize)

################################
# Variables
################################

plots_dir = here::here("plots", "introgression")

# HdrR chromosome lengths

chroms = read.table(here::here("data/Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt")) %>% 
  dplyr::select(chr = V1, end = V2) %>% 
  dplyr::filter(chr != "MT") %>% 
  dplyr::mutate(chr = paste("chr", chr, sep = ""),
                start = 0,
                end = as.numeric(end)) %>% 
  dplyr::select(chr, start, end)

################################
# Plotting
################################

# Factors

chr_order = c(seq(1,24), "all")
fish_order <- c("MIKK", "HdrR", "HSOK", "HNI")

# Palettes

pal_abba <- c("#9E2B25", "#F3B61F", "#631E68", "#F6673A", "#F33A56", "#55B6B0", "#621B00")
names(pal_abba) <- c("iCab", "HdrR", "HSOK", "HNI", "melastigma", "javanicus", "Kaga")
