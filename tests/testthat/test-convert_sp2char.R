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

test_that("convert_sp2char function outputs correct structure", {
    transform_type <- c('none', 'log10', 'square', 'cube')
    for (i in seq(transform_type)) {
        resultSE <- convert_sp2char(processed_se, transform=transform_type[i])
        expect_s4_class(resultSE, "SummarizedExperiment")
    }
})

test_that("convert_sp2char function outputs correct structure when rowData not contains '|'.", {
    rowData <- as.data.frame(SummarizedExperiment::rowData(processed_se))
    rowData <- rowData %>%
        dplyr::mutate(across(everything(), ~ifelse(grepl("\\|", .), NA, .)))
    SE <- processed_se
    SummarizedExperiment::rowData(SE) <- rowData
    resultSE <- convert_sp2char(processed_se, transform='log10')
    expect_s4_class(resultSE, "SummarizedExperiment")
})
