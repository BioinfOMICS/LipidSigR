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
    normalization='Percentage', transform='log10')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

## Multi Group
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_multiGroup <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')

#### Function ####
ora_test <- function(res, nGroup, significant, p_cutoff){
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_barPlot, "ggplot")
    expect_s3_class(res$interactive_barPlot, "plotly")
    # Verify if it is a data frame.
    expect_s3_class(res$enrich_result, 'data.frame')
    expect_s3_class(res$table_barPlot, 'data.frame')
    # Confirm the existence of required columns.
    expect_true(all(c('significance', 'aspect', 'characteristic', 'charFeature',
                      'direction', 'p', 'negLog10p', 'yText', 'yOrder', 'hover')
                    %in% colnames(res$table_barPlot)))
    if(nGroup == 'two'){
        expect_true(all(c('significance', 'aspect', 'characteristic', 'charFeature',
                          'condition', 'pval', 'padj', 'negLog10pval', 'negLog10padj',
                          'InCharDeLipid', 'notInCharDeLipid', 'InCharLipid',
                          'notInCharLipid', 'Lipids')
                        %in% colnames(res$enrich_result)))
    }else if(nGroup == 'multiple'){
        expect_true(all(c('significance', 'aspect', 'characteristic', 'charFeature',
                          'pval', 'padj', 'negLog10pval', 'negLog10padj',
                          'InCharDeLipid', 'notInCharDeLipid', 'InCharLipid',
                          'notInCharLipid', 'Lipids')
                        %in% colnames(res$enrich_result)))
    }
    # Confirm the correctness of significance.
    expect_equal(res$table_barPlot$charFeature[which(res$table_barPlot$p < p_cutoff)],
                 res$table_barPlot$charFeature[which(res$table_barPlot$significance == 'Yes')])
    # Confirm the correctness of significance.
    if(significant == 'pval'){
        expect_equal(res$enrich_result$charFeature[which(as.numeric(res$enrich_result$pval) < p_cutoff)],
                     res$enrich_result$charFeature[which(res$enrich_result$significance == 'Yes')])
    }else if(significant == 'padj'){
        expect_equal(res$enrich_result$charFeature[which(as.numeric(res$enrich_result$padj) < p_cutoff)],
                     res$enrich_result$charFeature[which(res$enrich_result$significance == 'Yes')])
    }
}

#### Test ####
test_that("Test under normal conditions: two groups", {
    for(significant in c('pval', 'padj')){
        res.all <- enrichment_ora(deSp_se=deSp_se_twoGroup, char=NULL,
                                  significant, p_cutoff)
        ora_test(res=res.all, nGroup='two', significant, p_cutoff)
        res.one <- enrichment_ora(deSp_se=deSp_se_twoGroup, char='class',
                                  significant, p_cutoff)
        ora_test(res=res.one, nGroup='two', significant, p_cutoff)
    }
    deSp_se <- deSp_twoGroup(
        processed_se_twoGroup, ref_group='ctrl', test='t-test',
        significant='padj', p_cutoff=0.05, FC_cutoff=1, transform='log10')
    res.all <- enrichment_ora(deSp_se, char=NULL, significant='pval', p_cutoff=0.05)
    ora_test(res=res.all, nGroup='two', significant='pval', p_cutoff=0.05)
})

test_that("Test under normal conditions: multiple groups", {
    for(significant in c('pval', 'padj')){
        res.all <- enrichment_ora(deSp_se=deSp_se_multiGroup, char=NULL,
                                  significant, p_cutoff)
        ora_test(res=res.all, nGroup='multiple', significant, p_cutoff)
        res.one <- enrichment_ora(deSp_se=deSp_se_multiGroup, char='class',
                                  significant, p_cutoff)
        ora_test(res=res.one, nGroup='multiple', significant, p_cutoff)
    }
})

test_that("Test for cases where there are no significant lipids.", {
    S4Vectors::metadata(deSp_se_multiGroup)$all_deSp_result <-
        S4Vectors::metadata(deSp_se_multiGroup)$all_deSp_result %>%
        dplyr::filter(pval > 0.1)
    expect_warning(
        enrichment_ora(
            deSp_se=deSp_se_multiGroup, char=NULL, significant='padj', p_cutoff),
        "Your data does not contain any significant lipids. Please reset the cutoff or review the differential expression results."
        )
})

test_that("Test for incorrect input and related error messages", {
    expect_error(
        enrichment_ora(deSp_se_twoGroup, char='lipid', significant='pval', p_cutoff=0.05),
        "The 'char' parameter must be one of the listed lipid characteristics. Please use the list_lipid_char() function to view which lipid characteristics can be analyzed."
    )
    expect_error(
        enrichment_ora(deSp_se_twoGroup, char='class', significant='Pval', p_cutoff=0.05),
        "The 'significant' parameter must be either 'pval' or 'padj'."
    )
    expect_error(
        enrichment_ora(deSp_se_twoGroup, char='class', significant='pval', p_cutoff=1.5),
        "The 'p_cutoff' parameter must be a numeric value between 0 and 1."
    )
})
