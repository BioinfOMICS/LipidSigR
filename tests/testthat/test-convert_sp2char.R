library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# load example data
data("de_data_twoGroup")
processed_se <- data_process(
     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
     replace_na_method='min', replace_na_method_ref=0.5,
     normalization='Percentage')

expect_char_abundance <- function(se, res){
    expect_s4_class(res, "SummarizedExperiment")
    expect_identical(SummarizedExperiment::colData(se),
                     SummarizedExperiment::colData(res))
    expect_false(identical(SummarizedExperiment::rowData(se),
                           SummarizedExperiment::rowData(res)))
}

test_that("convert_sp2char function outputs correct structure", {
    transform_type <- c('none', 'log10', 'square', 'cube')
    for (i in seq(transform_type)) {
        resultSE <- convert_sp2char(processed_se, transform=transform_type[i])
        expect_char_abundance(processed_se, resultSE)
    }
})

test_that("convert_sp2char function outputs correct structure when rowData not contains '|'.", {
    rowData <- as.data.frame(SummarizedExperiment::rowData(processed_se))
    rowData %<>%
        dplyr::mutate(across(everything(), ~ifelse(grepl("\\|", .), NA, .)))
    SE <- processed_se
    SummarizedExperiment::rowData(SE) <- rowData
    resultSE <- convert_sp2char(SE, transform='log10')
    expect_char_abundance(SE, resultSE)
})

test_that("convert_sp2char return error message for special condition or incorrect input", {
    expect_error(
        convert_sp2char(processed_se, transform=NULL),
        "transform must be one of 'none', 'log10', 'square' or 'cube'.")
    sub <- processed_se[1:2, ]
    resultSE <- convert_sp2char(sub, transform="log10")
    expect_char_abundance(sub, resultSE)
    rowData_sub <- as.data.frame(SummarizedExperiment::rowData(sub))
    rowData_sub[, colnames(rowData_sub[,-1])] <- NA
    SummarizedExperiment::rowData(sub) <- rowData_sub
    expect_error(convert_sp2char(sub, transform='log10'),
                 "No enough data.")
})
