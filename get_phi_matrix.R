#!/usr/local/bin/Rscript

options(scipen = 999)

source('kingraph.R')
file <- args <- commandArgs(TRUE)[1]

phi <- read.phi(file)

# Writing output to stdout
write.table(file = "", adj.from.phi(phi)$matrix, quote = FALSE, sep = "\t")

