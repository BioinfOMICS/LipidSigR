library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)

data("de_data_twoGroup")
data("se_multiGroup")
p_cutoff <- 0.05

## Two Group
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
## Multi Group
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_multiGroup <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')

#### Function ####
lsea_test <- function(res, significant, p_cutoff){
    # Verify if it is a list.
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_barPlot, "ggplot")
    expect_s3_class(res$interactive_barPlot, "plotly")
    # Verify if it is a data frame/list/double.
    expect_s3_class(res$enrich_result, 'data.frame')
    expect_s3_class(res$table_barPlot, 'data.frame')
    expect_type(res$lipid_set, 'list')
    expect_type(res$ranked_list, 'double')
    # Confirm the existence of required columns.
    expect_true(all(c('significance', 'aspect', 'characteristic', 'charFeature',
                      'p', 'NES', 'yText', 'hover')
                    %in% colnames(res$table_barPlot)))
    expect_true(all(c('significance', 'aspect', 'characteristic', 'charFeature',
                      'pval', 'padj', 'log2err', 'ES', 'NES', 'size', 'leadingEdge')
                    %in% colnames(res$enrich_result)))
    # Confirm the correctness of significance.
    if(significant == 'pval'){
        expect_equal(res$enrich_result$charFeature[which(res$enrich_result$pval < p_cutoff)],
                     res$enrich_result$charFeature[which(res$enrich_result$significance == 'Yes')])
    }else if(significant == 'padj'){
        expect_equal(res$enrich_result$charFeature[which(res$enrich_result$padj < p_cutoff)],
                     res$enrich_result$charFeature[which(res$enrich_result$significance == 'Yes')])
    }
}

#### Test ####
test_that("Test under normal conditions: two groups", {
    for(rank_by in c('log2FC', 'pval', 'padj', 'statistic')){
        for(significant in c('pval', 'padj')){
            res.all <- suppressWarnings(
                enrichment_lsea(
                    deSp_se=deSp_se_twoGroup, char=NULL, rank_by,
                    significant, p_cutoff, n_lipid=2))
            lsea_test(res=res.all, significant, p_cutoff)
            res.one <- suppressWarnings(
                enrichment_lsea(
                    deSp_se=deSp_se_twoGroup, char='class', rank_by,
                    significant, p_cutoff, n_lipid=2))
            lsea_test(res=res.one, significant, p_cutoff)
        }
    }
})

test_that("Test under normal conditions: multiple groups", {
    for(rank_by in c('pval', 'padj', 'statistic')){
        for(significant in c('pval', 'padj')){
            res.all <- suppressWarnings(
                enrichment_lsea(
                    deSp_se=deSp_se_multiGroup, char=NULL, rank_by,
                    significant, p_cutoff, n_lipid=2))
            lsea_test(res.all, significant, p_cutoff)
            res.one <- suppressWarnings(
                    enrichment_lsea(
                    deSp_se=deSp_se_multiGroup, char='class', rank_by,
                    significant, p_cutoff, n_lipid=2))
            lsea_test(res=res.one, significant, p_cutoff)
        }
    }
    expect_error(enrichment_lsea(
        deSp_se=deSp_se_multiGroup, char=NULL, rank_by='log2FC',
        significant='pval', p_cutoff, n_lipid=2),
        info='Log2 fold change (log2FC) ranking is not suitable for multiple group data. Please choose p-value (pval), adjusted p-value (padj), or statistic as the ranking criterion instead.')
})
