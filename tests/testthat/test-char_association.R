library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)

data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

p_cutoff <- 0.05
FC_cutoff <- 1


#### Function ####
char_association_test <- function(res, nGroup){
    ## Lollipop plot & Word cloud
    # Verify if it is a data frame.
    expect_true(is.data.frame(res$table_lollipop))
    expect_true(is.data.frame(res$table_wordCloud))
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_lollipop, "ggplot")
    expect_s3_class(res$interactive_lollipop, "plotly")
    expect_s3_class(res$static_wordCloud, "recordedplot")
    expect_s3_class(res$interactive_wordCloud, "htmlwidget")
    # Confirm the existence of required columns.
    expect_true(all(c('characteristic', 'freqs') %in% colnames(res$table_wordCloud)))
    if(nGroup == 'two'){
        # Confirm the existence of required columns.
        expect_true(all(c(
            'feature', 'FC', 'log2FC', 'characteristic', 'pval', 'negLog10pval',
            'padj', 'negLog10padj', 'hover', 'x') %in% colnames(res$table_lollipop)))
        ## Bar plot
        # Verify if it is a data frame.
        expect_true(is.data.frame(res$table_barPlot))
        # Explicitly force the evaluation of the plot.
        expect_s3_class(res$static_barPlot, "ggplot")
        expect_s3_class(res$interactive_barPlot, "plotly")
        # Confirm the existence of required columns.
        expect_true(all(c(
            'characteristic', 'log2FC.mean', 'log2FC.sd', 'log2FC.direction',
            'significant', 'hover', 'x') %in% colnames(res$table_barPlot)))
        # Confirm the correctness of significance.
        expect_equal(res$table_barPlot$characteristic[which(res$table_barPlot$significant == 'Yes')],
                     res$table_barPlot$characteristic[which(abs(res$table_barPlot$log2FC.mean) > log2(FC_cutoff))])
    }else if(nGroup == 'multiple'){
        # Confirm the existence of required columns.
        expect_true(all(c(
            'feature', 'characteristic', 'pval', 'negLog10pval', 'padj',
            'negLog10padj', 'y', 'hover', 'x') %in% colnames(res$table_lollipop)))
        expect_false(all(c('FC', 'log2FC') %in% colnames(res$table_lollipop)))
        ## Bar plot
        # Verify if it is a data frame.
        expect_equal(res$table_barPlot, "Only provided when two-group differential analysis data is input.")
        # Explicitly force the evaluation of the plot.
        expect_equal(res$static_barPlot, "Only provided when two-group differential analysis data is input.")
        expect_equal(res$interactive_barPlot, "Only provided when two-group differential analysis data is input.")
    }
}


#### Test ####
test_that("Test under normal conditions: two groups", {
    for(test in c('t-test', 'Wilcoxon test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                deSp_se_twoGroup <- deSp_twoGroup(
                    processed_se_twoGroup, ref_group='ctrl', test, significant,
                    p_cutoff, FC_cutoff, transform)
                res <- char_association(deSp_se_twoGroup, char='class')
                char_association_test(res, nGroup='two')
            }
        }
    }
})

test_that("Test under normal conditions: multi groups", {
    for(test in c('One-way ANOVA', 'Kruskal-Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                deSp_se_multiGroup <- deSp_multiGroup(
                    processed_se_multiGroup, ref_group='ctrl', test,
                    significant, p_cutoff, transform)
                res <- char_association(deSp_se_multiGroup, char='class')
                char_association_test(res, nGroup='multiple')
            }
        }
    }
})

test_that("Test for error char input", {
    deSp_se_multiGroup <- deSp_multiGroup(
        processed_se_multiGroup, ref_group='ctrl', test="One-way ANOVA",
        significant="pval", p_cutoff=0.005, transform="log10")
    expect_error(
        char_association(deSp_se_multiGroup, char=NULL),
        "Wrong char input, you can view the available char list by list_lipid_char function.")
})
