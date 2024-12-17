library(testthat)
# library(SummarizedExperiment)
# library(tidyverse)

data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

#### Function ####
ratioChar_se_format <- function(se) {
    if (!is.null(se)) {
        # Check if `se` is an SE object
        expect_true(inherits(se, "SummarizedExperiment"),
                    info = "Object is not an SummarizedExperiment")
        # Check if mandatory components are present
        expect_true("assays" %in% slotNames(se),
                    info = "Assays slot is missing")
        # expect_true("rowData" %in% slotNames(se),
        #             info = "rowData slot is missing")
        expect_true("colData" %in% slotNames(se),
                    info = "colData slot is missing")
        expect_true("metadata" %in% slotNames(se),
                    info = "metadata slot is missing")
        # Check the structure of assays, rowData, and colData
        # Check expression data exists and is a matrix
        expect_true(is.matrix(SummarizedExperiment::assay(se)),
                    info = "Expression data is not a matrix")
        expect_true(all(dim(SummarizedExperiment::assay(se)) > 0),
                    info = "Expression data has non-positive dimensions")
    }
}


#### Test ####
test_that("Test under normal conditions", {
    for(transform in c('none', 'log2')){
        ratio.res2 <- convert_sp2ratio(processed_se_twoGroup, transform)
        ratio.res3 <- convert_sp2ratio(processed_se_multiGroup, transform)
        ratioChar_se_format(se=ratio.res2)
        ratioChar_se_format(se=ratio.res3)
    }
})


test_that("Test without any molecular subspecies", {
    ## Without any molecular subspecies
    processed_se_multiGroup <- processed_se_multiGroup[19:34,]
    for(transform in c('none', 'log2')){
        ratio.res3 <- suppressWarnings(
            convert_sp2ratio(processed_se_multiGroup, transform))
        expect_true(is.null(ratio.res3))
        expect_warning(convert_sp2ratio(processed_se_multiGroup, transform),
                       info='There are no ratio characteristics that can be converted in your dataset.')
    }
})

test_that("Test lipids containing only even-chain fatty acids", {
    ## Lipids containing only even-chain fatty acids.
    processed_se_multiGroup <- processed_se_multiGroup[72:202,]
    processed_se_multiGroup <- processed_se_multiGroup[c(1:2, 5:10, 19:21, 70, 72:73, 75:86),]
    for(transform in c('none', 'log2')){
        ratio.res3 <- convert_sp2ratio(processed_se_multiGroup, transform)
        ratioChar_se_format(se=ratio.res3)
        expect_message(convert_sp2ratio(processed_se_multiGroup, transform),
                       info='There are 1 ratio characteristics that can be converted in your dataset.')
    }
})

test_that("Test for special condition.", {
    sub <- processed_se_twoGroup[1:2, ]
    expect_warning(
        convert_sp2ratio(sub, transform="log2"),
        "There are no ratio characteristics that can be converted in your dataset.")
    res <- suppressWarnings(
        convert_sp2ratio(sub, transform="log2")
    )
    ratioChar_se_format(res)
})

test_that("Test for error input", {
    expect_error(
        convert_sp2ratio(processed_se_multiGroup, transform="log10"),
        "The 'transform' parameter must be either 'none' or 'log2'.")
})
