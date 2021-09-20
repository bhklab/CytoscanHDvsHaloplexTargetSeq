# 0 -- Load dependencies
renv::activate()

library(GenomicRanges)
library(data.table)
library(qs)
library(wCorr) # for weighted correlations
library(SNFtool) # for calNMI
library(aricode) # for NMI
library(Metrics) # for RMSE
library(ggplot2)


# 0 -- Parse Snakemake arguments
# input <- snakemake@input
input <- list(
    cytoscan=list.files('results', pattern='cyt.*grlist.qs', full.names=TRUE),
    haloplex=list.files('results', pattern='halo.*_segment_grlist.qs', 
        full.names=TRUE)
)
# output <- snakemake@output
output <- 'results'

## 1 ---- Read in the GRangesList objects
grlist_list <- lapply(input, qread)

cytoscan_grlist <- grlist_list$cytoscan
haloplex_grlist <- grlist_list$haloplex
names(haloplex_grlist) <- 
    gsub('haloplex|_|.cnr|.cns', '', names(haloplex_grlist))

## 1.1 -- Load some metadata to help match samples
params <- list(cytoscan_meta='metadata/fletcher_LMS_pairs.tsv')
cytoscan_meta <- fread(params$cytoscan_meta)
cel_file_names <- gsub('^rawdata\\/|.CEL$', '', cytoscan_meta$CEL)
names(cytoscan_grlist) <- cel_file_names

# 1.2 -- Find matching sample between the array and targets seq data
# There is a case where it looks like LMS90 was merged with LMS92? 
# I will add this manually for now.
haloplex_grlist[['LMS_90-422_92-543']] <- c(haloplex_grlist[['LMS90-422']], 
    haloplex_grlist[['LMS92-543']])
shared_samples <- intersect(names(cytoscan_grlist), names(haloplex_grlist))

## 2 ---- Merge by overlapping segments

olaps_list <- as(mapply(nearest, x=haloplex_grlist[shared_samples], 
    subject=cytoscan_grlist[shared_samples]), 'List')
cytoscan_olaps_segment_list <- structure(mendoapply(`[`, 
    x=cytoscan_grlist[shared_samples], i=olaps_list), .Names=shared_samples)
haloplex_olaps_segment_list <- haloplex_grlist[shared_samples]

# merge metadata
haloplex_flat <- unlist(haloplex_olaps_segment_list)
cytoscan_flat <- unlist(cytoscan_olaps_segment_list)
mcols(haloplex_flat)$cytoscan_log2 <- cytoscan_flat$seg.mean
mcols(haloplex_flat)$cytoscan_cn <- cytoscan_flat$TCN
mcols(cytoscan_flat)$haloplex_seg.mean <- haloplex_flat$log2

# relist
haloplex_olaps_segment_list <- relist(haloplex_flat, haloplex_olaps_segment_list)
cytoscan_olaps_segment_list <- relist(cytoscan_flat, cytoscan_olaps_segment_list)

# convert to data.tables
haloplex_dt <- as.data.table(haloplex_olaps_segment_list)
cytoscan_dt <- as.data.table(cytoscan_olaps_segment_list)

# drop sex chromosomes
haloplex_dt <- haloplex_dt[!(as.vector(seqnames) %like% 'chrY|chrX'), ]

# rename the sample column
setnames(haloplex_dt, c('seqnames', 'group_name'), c('chrom', 'sample'))

fwrite(haloplex_dt, file.path(output, 'haloplex_to_cytoscan_overlaps.csv'))

## 3 -- Calculate a similarity metrics based on 

# 3.1 -- Calculate weighted and unweighted similarity metrics

## FIXME: (A) Pearson correlation is not a valid statistic here due to:
##   (1) Violation of assumption of normality
##   (2) Presence of outliers
##   (...) There are probaly more assumptions violated

try_cor.test_pval <- function(...) tryCatch({ cor.test(...)$p.value }, 
    error=function(e) NA_real_)

# comparison data.table
haloplex_comparison <- haloplex_dt[,
    .(
        log2=mean(log2, na.rm=TRUE),
        pct_match_tcn=sum(cn == cytoscan_cn)/length(cn),
        ratio_diff_lt0.5=sum(abs(log2 - cytoscan_log2) < 0.5)/length(log2),
        euc_dist=dist(rbind(log2, cytoscan_log2)),
        spearman_cor=cor(log2, cytoscan_log2, method='spearman'),
        spearman_wcor=weightedCorr(log2, cytoscan_log2, weights=probes, 
            method='Spearman'),
        spearman_pvalue=try_cor.test_pval(log2, cytoscan_log2, 
            method='spearman'),
        nmi_tcn=NMI(as.character(cn), as.character(cytoscan_cn)),
        log2_over_rmse=mean(abs(log2/rmse(log2, cytoscan_log2)), na.rm=TRUE),
        read_depth=mean(depth, na.rm=TRUE),
        num_probes=mean(probes, na.rm=TRUE)
        num_segments=.N
    ),
    by=.(sample, chrom)]

# Drop NA correlations
haloplex_comparison <- haloplex_comparison[
    !is.na(spearman_cor), ]

haloplex_comparison <- rbindlist(
    list(haloplex_comparison,
        haloplex_dt[,
            .(
                chrom='chrAll',
                log2=mean(log2, na.rm=TRUE),
                pct_match_tcn=sum(cn == cytoscan_cn)/length(cn),
                ratio_diff_lt0.5=sum(abs(log2 - cytoscan_log2) < 0.5)/length(log2),
                euc_dist=dist(rbind(log2, cytoscan_log2)),
                spearman_cor=cor(log2, cytoscan_log2, method='spearman'),
                spearman_wcor=weightedCorr(log2, cytoscan_log2, weights=probes, 
                    method='Spearman'),
                spearman_pvalue=try_cor.test_pval(log2, cytoscan_log2, 
                    method='spearman'),
                nmi_tcn=NMI(as.character(cn), as.character(cytoscan_cn)),
                log2_over_rmse=mean(abs(log2/rmse(log2, cytoscan_log2)), na.rm=TRUE),
                read_depth=mean(depth, na.rm=TRUE),
                num_probes=mean(probes, na.rm=TRUE),
                num_segments=.N
            ),
            by=sample]),
    use.names=TRUE)

# save the results
fwrite(haloplex_comparison, file.path(output, 
    'haloplex_to_cytoscan_comparison.csv'))

# NOTE: Removed, see FIXME (A)
# -- 3.2 -- Visualize the similarity comparisons

heatmap_metric <- function(column, midpoint) {

    plot <- ggplot(haloplex_comparison, aes(x=sample, y=chrom)) +
        geom_tile(aes_string(fill=column)) +
        scale_fill_gradient2(low='blue', mid='white', high='red') +
        geom_text(aes_string(label=paste0('round(', column, ', 2)'))) +
        ggtitle(column)
    
    if (column %in% c('spearman_pvalue', )) {
        plot <- plot + scale_fill_gradient2(low='red', mid='white', high='blue')
    }
    
    return(plot)
}

heatmap_list <- mapply(FUN=heatmap_metric,
    column=colnames(haloplex_comparison[, .SD[, -'log2'], .SDcol=is.numeric]),
    midpoint=c(0.5, 0, 0, 0, 0.5, 1),
    SIMPLIFY=FALSE
    )

pdf(file.path(output, 'figures', 'haloplex_to_cytoscan_heatmaps.pdf'))
for (plot in heatmap_list) {
    print(plot)
}
dev.off()

