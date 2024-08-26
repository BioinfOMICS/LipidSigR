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

#### Test ####
test_that("Test under normal conditions: two groups", {
    lsea_res <- enrichment_lsea(
        deSp_se=deSp_se_twoGroup, char='class', rank_by='log2FC',
        significant='pval', p_cutoff, n_lipid=2)
    res <- plot_enrichment_lsea(lsea_res, char='class', char_feature='PC')
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_enrichPlot, "ggplot")
    expect_s3_class(res$interactive_enrichPlot, "plotly")
})

test_that("Test under normal conditions: multiple groups", {
    lsea_res <- suppressWarnings(
        enrichment_lsea(
            deSp_se=deSp_se_multiGroup, char='class', rank_by='pval',
            significant='pval', p_cutoff, n_lipid=2))
    res <- plot_enrichment_lsea(lsea_res, char='class', char_feature='PC')
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_enrichPlot, "ggplot")
    expect_s3_class(res$interactive_enrichPlot, "plotly")
})
