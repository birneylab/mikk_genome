#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)
library(ape)

# Collect variables from bash script

file_in <- args[1]
dir_out <- args[2]

# Paste into new names

file_data <- paste(file_in, ".data.txt", sep = "")
file_tree <- paste(file_in, ".tree.txt", sep = "")
file_seq <- paste(file_in, ".seq.txt", sep = "")

# Get tree

phylo_tree <- ape::read.tree(file = file_tree)

# Find most recent common ancestor

ids <- phylo_tree$tip.label[grep("Olat|Ohni|Ohso|Omel|Ojav", phylo_tree$tip.label)] # get vector of IDs
node_number <- ape::getMRCA(phylo_tree, tip = ids) # get node number of most recent common ancestor
mrca_label <- c(phylo_tree$tip.label, phylo_tree$node.label)[node_number] # get label of node
mrca_label_trim <- stringr::str_split(mrca_label, "_", simplify = T)[2:4] %>% str_c(collapse = "_") # trim to match line in SEQ file

# Read in SEQ data

seq_dat <- readr::read_delim(file_seq,
                             delim = " ",
                             col_names = F)

# Re-jig names so they're in the same vector and get names and cols for the DATA file

org_names <- ifelse(grepl("ancestral", seq_dat$X2), seq_dat$X3, seq_dat$X2)
target_indices <- sort(c(grep(mrca_label_trim, org_names), grep("oryzias", org_names))) # get indices
target_names <- org_names[target_indices]

# Read in DATA file

dat <- data.table::fread(file_data,
                         header = F,
                         select = target_indices) # filter for target columns
colnames(dat) <- target_names # set column names
dat <- dat[dat$oryzias_latipes != "-", ] # remove all rows with "-" in HdrR column

# Get sequence start and end coords

chr <- stringr::str_split(basename(file_in), pattern = "_", simplify = T)[1]
seq_start <- stringr::str_split(basename(file_in), pattern = "_", simplify = T)[2]
seq_end <- stringr::str_split(basename(file_in), pattern = "_", simplify = T)[3]

# Add to data frame and re-order

dat$chr <- chr
dat$coord <- seq(seq_start, seq_end)
dat <- dat %>% dplyr::select(chr, coord, everything())

# Save to file

bname_out <- paste(basename(file_in), ".txt", sep = "")
file_out <- file.path(dir_out, bname_out)
write.table(dat, file_out, quote = F, sep = "\t", row.names = F)
