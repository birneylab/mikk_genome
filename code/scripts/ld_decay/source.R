#!/usr/bin/env Rscript

library(here)
library(tidyverse)
library(gaston)
library(heatmaply)
library(biomaRt)
library(ape)
library(cowplot)
library(plotly)

# To avoid the following biomaRt error:
#Ensembl site unresponsive, trying asia mirror
#Error in curl::curl_fetch_memory(url, handle = handle) : 
#  SSL certificate problem: unable to get local issuer certificate
httr::set_config(httr::config(ssl_verifypeer = FALSE))

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

# Big data storage on Codon
lts_dir = "/nfs/research/birney/users/ian/mikk_genome/ld_decay"
