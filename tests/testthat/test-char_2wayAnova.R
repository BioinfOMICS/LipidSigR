library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)

data("de_data_twoGroup")
data("se_multiGroup")
processed_se_twoGroup <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

#### Function ####
twoWayAnova_test <- function(res){
    # Verify if it is a data frame.
    expect_true(is.data.frame(res))
    # Confirm whether columns that cannot be NA indeed do not contain any NA values.
    expect_true(all(!is.na(res$aspect)))
    expect_true(all(!is.na(res$characteristic)))
    # Confirm whether the data frame consists of numeric values.
    expect_true(all(apply(res[,-1:-2], 2, function(x) all(is.numeric(x)))))
    # Confirm the existence of required columns.
    expect_true(all(c(
        "aspect", "characteristic", "fval_2factors", "pval_2factors",
        "padj_2factors", "fval_feature", "pval_feature", "padj_feature",
        "fval_group", "pval_group", "padj_group") %in% colnames(res)))
}


#### Test ####
test_that("Test under normal conditions: 2 groups", {
    for(ratio_transform in c('none', 'log2')){
        for(char_transform in c('none', 'log10', 'square', 'cube')){
            res <- char_2wayAnova(processed_se=processed_se_twoGroup,
                                  ratio_transform, char_transform)
            twoWayAnova_test(res)
        }
    }
})


test_that("Test under normal conditions: 3 groups", {
    for(ratio_transform in c('none', 'log2')){
        for(char_transform in c('none', 'log10', 'square', 'cube')){
            res <- char_2wayAnova(processed_se=processed_se_multiGroup,
                                  ratio_transform, char_transform)
            twoWayAnova_test(res)
        }
    }
})




