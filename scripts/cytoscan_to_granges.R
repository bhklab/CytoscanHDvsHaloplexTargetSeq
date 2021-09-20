# 0 -- Load dependencies
renv::activate()

library(GenomicRanges)
library(qs)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

input <- list.files('procdata', '*grList.qs', full.names=TRUE)
output <- 'results'

granges_list <- qread(input)

qsave(granges_list, file.path(output, 'cytoscan_grlist.qs'))