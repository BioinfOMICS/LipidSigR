library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("abundance_twoGroup")
data("group_info_twoGroup")
two_parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_twoGroup$feature)
two_recognized_lipid <- two_parse_lipid$Original.Name[which(two_parse_lipid$Grammar != 'NOT_PARSEABLE')]
two_abundance <- abundance_twoGroup %>%
    dplyr::filter(feature %in% two_recognized_lipid)
two_char_table <- two_parse_lipid %>%
    dplyr::filter(Original.Name %in% two_recognized_lipid)
se <- as_summarized_experiment(
    two_abundance, two_char_table, group_info_twoGroup, n_group="two",
    paired_sample=FALSE)

se_error <- SummarizedExperiment::SummarizedExperiment(
    assays=NULL, rowData=NULL, colData=NULL)

expect_correct_se <- function(
        input_se, processed_se, exclude_missing, replace_na_method, normalization) {
    expect_true(inherits(processed_se, "SummarizedExperiment"))
    expect_true("assays" %in% slotNames(processed_se),
                info = "Assays slot is missing")
    expect_true("colData" %in% slotNames(processed_se),
                info = "colData slot is missing")
    expect_true("elementMetadata" %in% slotNames(processed_se),
                info = "elementMetadata slot is missing")
    ## expect value
    # expect_identical(
    #     SummarizedExperiment::colData(input_se), SummarizedExperiment::colData(processed_se))
    # if (exclude_missing==FALSE && replace_na_method=='none' && normalization=='none') {
    #     expect_identical(
    #         SummarizedExperiment::assay(input_se), SummarizedExperiment::assay(processed_se))
    # } else {
    #     expect_false(
    #         identical(SummarizedExperiment::assay(input_se), SummarizedExperiment::assay(processed_se)))
    # }
}


# Test for basic functionality and output structure
test_that("data_process function work correctly", {
    processed_se <- data_process(
        se, exclude_missing=TRUE, exclude_missing_pct=70,
        replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage')
    expect_correct_se(se, processed_se, exclude_missing=TRUE, replace_na_method='min', normalization='Percentage')
    expect_error(data_process(
        se, exclude_missing=FALSE, exclude_missing_pct=70,
        replace_na_method='none', replace_na_method_ref=0.5, normalization='none'),
        "Detect 0 values or NAs in abundance data. Please select a data imputation method in replace_na_method other than 'none'.")
    expect_error(data_process(
        se_error, exclude_missing=TRUE, exclude_missing_pct=70, replace_na_method='min',
        replace_na_method_ref=0.5, normalization='Percentage'))
})


## if all been removed....
