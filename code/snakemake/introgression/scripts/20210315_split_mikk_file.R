#!/usr/bin/env Rscript

samples <- scan(snakemake@input[[1]], what = "chr")
set.seed(54)
pop_1 <- sample(samples, 32)
pop_2 <- samples[!samples %in% pop_1]

write(samples, snakemake@output[[1]])
write(pop_1, snakemake@output[[2]])
write(pop_2, snakemake@output[[3]])
