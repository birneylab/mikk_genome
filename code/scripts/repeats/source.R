#!/usr/bin/env Rscript

#############################
# Libraries
#############################

library(here)
library(tidyverse)
library(cowplot)
library(scales)
library(karyoploteR)
library(circlize)

#############################
# Functions
#############################

round.choose <- function(x, roundTo, dir = 1) {
  if(dir == 1) {  ##ROUND UP
    x + (roundTo - x %% roundTo)
  } else {
    if(dir == 0) {  ##ROUND DOWN
      x - (x %% roundTo)
    }
  }
}

#############################
# Plotting
#############################

## Factor order

chr_order = paste("chr", 1:24, sep = "")
#class_order = c("Simple_repeat", "SINE", "LTR", "Satellite", "DNA", "LINE", "Unknown", "RC", "tRNA", "rRNA", "misc", "Retroposon", "snRNA")
strand_order = c("+", "-")

## Palettes

pal_electroangler = c("#7400b8", "#6930c3", "#5e60ce", "#5390d9", "#4ea8de", "#48bfe3", "#56cfe1", "#64dfdf", "#72efdd", "#80ffdb")
