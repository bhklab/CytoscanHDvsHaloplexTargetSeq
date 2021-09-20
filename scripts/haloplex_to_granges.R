# 0 -- Load dependencies
renv::activate()

library(GenomicRanges)
library(data.table)
library(qs)

# 0 -- Parse Snakemake arguments
# input <- snakemake@input
# params <- snakemake@params
# nthreads <- snakemake@threads
# output <- snakemake@output

# 1 -- Configure file paths
input <- list(
    segment=list.files('procdata', pattern='haloplex.*call.cns', full.names=TRUE),
    region=list.files('procdata', pattern='haloplex.*cnr', full.names=TRUE)
)
output <- 'results'

# 2 -- Load in segment and region data

# segment
segment_dt_list <- lapply(input$segment, fread)
granges_segment_list <- as(
    lapply(segment_dt_list, makeGRangesFromDataFrame, keep.extra.columns=TRUE), 
    'GRangesList')
names(granges_segment_list) <- vapply(input$segment, basename, character(1))
names(granges_segment_list) <- gsub('haloplex|_|.call.cns', '', 
    names(granges_segment_list))

qsave(granges_segment_list, file.path(output, 'haloplex_segment_grlist.qs'))

# region
region_dt_list <- lapply(input$region, fread)
granges_region_list <- as(
    lapply(region_dt_list, makeGRangesFromDataFrame, keep.extra.columns=TRUE), 
    'GRangesList')
names(granges_region_list) <- vapply(input$region, basename, character(1))
names(granges_region_list) <- gsub('haloplex|_|.cnr', '', 
    names(granges_region_list))

qsave(granges_region_list, file.path(output, 'haloplex_region_grlist.qs'))