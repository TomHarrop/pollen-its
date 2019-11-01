#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(dada2)

fwd_in <- snakemake@input[["r1"]]
rev_in <- snakemake@input[["r2"]]
fwd_out <- snakemake@output[["r1"]]
rev_out <- snakemake@output[["r2"]]

filterAndTrim(fwd_in,
              fwd_out,
              rev_in,
              rev_out,
              maxN = 0,
              multithread = FALSE,
              matchIDs = TRUE)

sessionInfo()
