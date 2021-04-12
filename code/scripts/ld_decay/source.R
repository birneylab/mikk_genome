#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(gaston)
library(heatmaply)
library(biomaRt)
library(ape)
library(cowplot)
library(plotly)


#####################
# Variables
#####################

chroms = read.table(here::here("data/Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt")) %>% 
  dplyr::select(chr = V1, end = V2) %>% 
  dplyr::filter(chr != "MT") %>% 
  dplyr::mutate(chr = paste("chr", chr, sep = ""),
                start = 0,
                end = as.numeric(end)) %>% 
  dplyr::select(chr, start, end)