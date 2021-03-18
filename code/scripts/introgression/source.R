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

################################
# Variables
################################

plots_dir = here::here("plots", "introgression")

################################
# Plotting
################################

# Factors

chr_order = c(seq(1,24), "all")
fish_order <- c("MIKK", "HdrR", "HSOK", "HNI")

# Palettes

pal_abba <- c("#F3B61F", "#631E68", "#F6673A", "#F33A56", "#55B6B0")
names(pal_abba) <- c("HdrR", "HSOK", "HNI", "melastigma", "javanicus")