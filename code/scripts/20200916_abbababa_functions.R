#!/usr/bin/Rscript

# D statistic
D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
    }

# F statistic
f.stat <- function(p1, p2, p3a, p3b) {
    ABBA_numerator <- (1 - p1) * p2 * p3a
    BABA_numerator <- p1 * (1 - p2) * p3a

    ABBA_denominator <- (1 - p1) * p3b * p3a
    BABA_denominator <- p1 * (1 - p3b) * p3a

    (sum(ABBA_numerator) - sum(BABA_numerator)) /
    (sum(ABBA_denominator) - sum(BABA_denominator))
    }

# Get block indices
get_block_indices <- function(block_size, positions, chromosomes = NULL){
    if (is.null(chromosomes) == TRUE) {
        block_starts <- seq(min(positions), max(positions), block_size)

        block_ends <- block_starts + block_size - 1

        lapply(1:length(block_starts), function(x) which(positions >= block_starts[x] &
                                                         positions <= block_ends[x]))
        }
    else {
        chrom_names <- unique(chromosomes)

        block_starts <- lapply(chrom_names, function(chrom_name) seq(min(positions[chromosomes==chrom_name]),
                                                                     max(positions[chromosomes==chrom_name]), block_size))

        block_chroms <- unlist(lapply(1:length(block_starts), function(x) rep(chrom_names[x], length(block_starts[[x]]))))

        block_starts <- unlist(block_starts)

        block_ends <- block_starts + block_size - 1

        lapply(1:length(block_starts), function(x) which(chromosomes == block_chroms[x] &
                                                     positions >= block_starts[x] &
                                                     positions <= block_ends[x]))
        }
    }

#this function runs the jackknife procedure by calculating pseudovalues by removing one block at a time
#if the arguments specified by "..." are vectors, they will be indexed as they are.
#if they have two dimensions, they will be indexed along the first dimension
get_jackknife_sd <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))))
        }
    else{
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))))
        }
    }

block_jackknife <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))
        }
    else{
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))
        }

    mean <- mean(pseudovalues)

    std_dev <- sd(pseudovalues)

    list(mean=mean, std_dev=std_dev)
    }

run_abbababa <- function(data, P1, P2){

  # set popns and remove NAs
  P3 <- "mikk"
  P3a <- "mikk_a"
  P3b <- "mikk_b"
  pops <- c(P1, P2, P3, P3a, P3b)

  # Select only those populations and remove NAs
  freq_table <- data %>%
    dplyr::select(chr, pos, ancestral, derived, all_of(pops)) %>%
    tidyr::drop_na()

  # Create output list
  out_list <- list()

  # Add populations
  out_list[["Populations"]] <- list("P1" = P1,
                                    "P2" = P2,
                                    "P3" = "mikk")

  # Get genome-wide D stat
  out_list[["Genome-wide D"]] <- list()

  D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])
  out_list[["Genome-wide D"]][["D statistic"]] <- D

  block_indices <- get_block_indices(block_size=1e6, # Block jackknife to obtain SD and Z-score
                                     positions=freq_table$pos,
                                     chromosomes=freq_table$chr)

  block_indices <- block_indices[lapply(block_indices, length) > 0 ]  # remove empty entries

  n_blocks <- length(block_indices)

  out_list[["Genome-wide D"]][["Number of blocks"]] <- n_blocks

  D_sd <- get_jackknife_sd(block_indices=block_indices, # get D SD
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])

  out_list[["Genome-wide D"]][["Standard deviation"]] <- D_sd

  D_err <- D_sd/sqrt(n_blocks)

  D_Z <- D / D_err

  out_list[["Genome-wide D"]][["Z score"]] <- D_Z

  # Get admixture proportion
  out_list[["Admixture"]] <- list()

  f <- f.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])

  out_list[["Admixture"]][["f statistic"]] <- f

  f_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=f.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])

  out_list[["Admixture"]][["Standard deviation"]] <- f_sd

  f_err <- f_sd/sqrt(n_blocks)

  f_CI_lower <- f - 1.96*f_err
  f_CI_upper <- f + 1.96*f_err

  out_list[["Admixture"]][["Confidence interval"]] <- list("lower" = f_CI_lower,
                                                           "upper" = f_CI_upper)

  # Get per-chromosome stats
  chrom_names <- unique(freq_table$chr)
  chrom_indices <- lapply(chrom_names, function(chrom) which(freq_table$chr == chrom))
  names(chrom_indices) <- chrom_names

  out_list[["Per-chromosome"]] <- lapply(chrom_names, function(chrom){
    per_chr_out <- list()
    # get D stat
    D_by_chrom <- D.stat(freq_table[chrom_indices[[chrom]], P1],
                         freq_table[chrom_indices[[chrom]], P2],
                         freq_table[chrom_indices[[chrom]], P3])

    per_chr_out[["D statistic"]] <- D_by_chrom

    # number of SNPs
    per_chr_out[["Number of SNPs"]] <- length(chrom_indices[[chrom]])

    # get block indices
    block_indices_by_chrom <- get_block_indices(block_size=1e6,
                                       positions=freq_table$pos[freq_table$chr == chrom])
    # remove empty blocks
    block_indices_by_chrom <- block_indices_by_chrom[lapply(block_indices_by_chrom, length) > 0 ]

    # Get SD
    D_sd_by_chrom <- get_jackknife_sd(block_indices=block_indices_by_chrom,
                                      FUN=D.stat,
                                      freq_table[chrom_indices[[chrom]], P1],
                                      freq_table[chrom_indices[[chrom]], P2],
                                      freq_table[chrom_indices[[chrom]], P3])

    per_chr_out[["Standard error"]] <- D_sd_by_chrom

    # Get Z score
    D_err_by_chrom <- D_sd_by_chrom / sqrt(length(block_indices_by_chrom))

    D_Z_by_chrom <- D_by_chrom / D_err_by_chrom

    per_chr_out[["Z score"]] <- D_Z_by_chrom

    # Get admixture
    per_chr_out[["Admixture"]] <- list()

    f_by_chrom <- f.stat(freq_table[chrom_indices[[chrom]],P1],
                         freq_table[chrom_indices[[chrom]],P2],
                         freq_table[chrom_indices[[chrom]],P3a],
                         freq_table[chrom_indices[[chrom]],P3b])

    per_chr_out[["Admixture"]][["f statistic"]] <- f_by_chrom

    f_sd_by_chrom <- get_jackknife_sd(block_indices=block_indices_by_chrom,
                                      FUN=f.stat,
                                      freq_table[chrom_indices[[chrom]],P1],
                                      freq_table[chrom_indices[[chrom]],P2],
                                      freq_table[chrom_indices[[chrom]],P3a],
                                      freq_table[chrom_indices[[chrom]],P3b])

    per_chr_out[["Admixture"]][["Standard deviation"]] <- f_sd_by_chrom

    f_err_by_chrom <- f_sd_by_chrom / sqrt(length(block_indices_by_chrom))

    f_CI_lower_by_chrom <- f - 1.96*f_err_by_chrom
    f_CI_upper_by_chrom <- f + 1.96*f_err_by_chrom

    per_chr_out[["Admixture"]][["Confidence interval"]] <- list("lower" = f_CI_lower_by_chrom,
                                                                "upper" = f_CI_upper_by_chrom)

    # Return list
    return(per_chr_out)
  })

  names(out_list[["Per-chromosome"]]) <- chrom_names

  # return output
  return(out_list)
}
