#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

# Import libraries

library(tidyverse)
library(ape)

# Collect variables from bash script

file_in <- args[1]
file_out <- args[2]

# Paste into new names

file_pref = stringr::str_remove(file_in, ".data.txt")
file_data <- paste(file_pref, ".data.txt", sep = "")
file_tree <- paste(file_pref, ".tree.txt", sep = "")
file_seq <- paste(file_pref, ".seq.txt", sep = "")

# Get sequence start and end coords

chr <- stringr::str_split(basename(file_pref), pattern = "_", simplify = T)[1]
seq_start <- stringr::str_split(basename(file_pref), pattern = "_", simplify = T)[2]
seq_end <- stringr::str_split(basename(file_pref), pattern = "_", simplify = T)[3]

# Get tree

phylo_tree <- ape::read.tree(file = file_tree)

# Find most recent common ancestor

ids <- phylo_tree$tip.label[grep("Olat|Ohni|Ohso", phylo_tree$tip.label)] # get vector of IDs
node_number <- ape::getMRCA(phylo_tree, tip = ids) # get node number of most recent common ancestor
mrca_label <- c(phylo_tree$tip.label, phylo_tree$node.label)[node_number] # get label of node
mrca_label_trim <- stringr::str_split(mrca_label, "_", simplify = T)[2:4] %>% str_c(collapse = "_") # trim to match line in SEQ file

# Read in SEQ data

seq_dat <- readr::read_delim(file_seq,
                             delim = " ",
                             col_names = F)

# Re-jig names so they're in the same vector and get names and cols for the DATA file

org_names <- ifelse(grepl("ancestral", seq_dat$X2), seq_dat$X3, seq_dat$X2)

# Get indices

# add target ancestor
keep_anc <- which(seq_dat$X3 %in% mrca_label_trim)
# add target HdrR
keep_hdrr <- which(seq_dat$X2 == "oryzias_latipes" & seq_dat$X4 == seq_start & seq_dat$X5 == seq_end)
# add target others
target_species <- c("oryzias_latipes_hni",
                    "oryzias_latipes_hsok")
# get their indices
other_inds <- unlist(sapply(target_species, function(x){
  target_ind <- which(seq_dat$X2 %in% x)
  # if there are multiple hits, take the one on the forward strand
  if (length(target_ind) > 1 ){
    target_ind <- which(seq_dat$X2 %in% x & seq_dat$X6 == "1")
    # if there are still multiple hits, take the first one
    if (length(target_ind) > 1){
      target_ind <- which(seq_dat$X2 %in% x & seq_dat$X6 == "1")[1]
    }
  }
  return(target_ind)
}))

target_indices <- c(keep_anc, keep_hdrr, other_inds)
target_indices <- unname(target_indices)

target_names <- org_names[target_indices]

# Read in DATA file

dat <- data.table::fread(file_data,
                         header = F,
                         select = target_indices) # filter for target columns
colnames(dat) <- target_names # set column names
dat <- dat[dat$oryzias_latipes != "-", ] # remove all rows with "-" in HdrR column

# Add to data frame and re-order

dat$chr <- chr
dat$coord <- seq(seq_start, seq_end)
dat <- dat %>% dplyr::select(chr,
                             coord,
                             oryzias_latipes,
                             everything())

# Save to file

write.table(dat, file_out, quote = F, sep = "\t", row.names = F)
