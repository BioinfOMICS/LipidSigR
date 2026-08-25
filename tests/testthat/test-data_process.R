library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("de_data_twoGroup")
se <- de_data_twoGroup

expect_correct_se <- function(
        input_se, processed_se, exclude_missing, replace_na_method, normalization, transform) {
    expect_s4_class(processed_se, "SummarizedExperiment")
    expect_true("assays" %in% slotNames(processed_se),
                info="Assays slot is missing")
    expect_true("colData" %in% slotNames(processed_se),
                info="colData slot is missing")
    expect_true("elementMetadata" %in% slotNames(processed_se),
                info="elementMetadata slot is missing")
    ## expect value
    expect_identical(
        SummarizedExperiment::colData(input_se), SummarizedExperiment::colData(processed_se))
    if (exclude_missing == FALSE && replace_na_method == 'none' && normalization == 'none' && transform =='none') {
        expect_identical(
            SummarizedExperiment::assay(input_se), SummarizedExperiment::assay(processed_se))
    }
}


# Test for basic functionality and output structure
test_that("data_process function work correctly", {
    processed_se <- data_process(
        se, exclude_missing=TRUE, exclude_missing_pct=70,
        replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
    expect_correct_se(se, processed_se, exclude_missing=TRUE, replace_na_method='min', normalization='Percentage', transform='log10')
    expect_snapshot(data_process(
        se, exclude_missing=FALSE, exclude_missing_pct=70,
        replace_na_method='none', replace_na_method_ref=0.5, normalization='none'), error=TRUE)
})

# Test for all combinations of imputation and normalization methods
test_that("data_process runs successfully on all combinations of imputation and normalization methods", {
    imputation_methods <- c(
        "BPCA", "GMS", "KNN", "impSeq","Impseqrob", "LLS", "mean", "median",
        "min", "MinProb", "PPCA", "PI", "QRILC", "RandomForest", "SVD")
    normalization_methods <- c('CubicSplines', 'CyclicLoess', 'EigenMS', 'RLR', 'nonLinear',
        'linear', 'Mean', 'Median', 'Percentage', 'PQN', 'Quantile', 'Sum', 'TMM' ,'none')
    transform_methods <- c('none', 'log10', 'cube', 'square')
    for (imp_method in imputation_methods) {
        for (norm_method in normalization_methods) {
            for (tran_methods in transform_methods) {
                processed_se <- suppressWarnings(suppressMessages(
                    data_process(
                        se, exclude_missing=TRUE, exclude_missing_pct=70,
                        replace_na_method=imp_method,
                        replace_na_method_ref=ifelse(
                            imp_method %in% c('QRILC', 'min'), 0.5,
                            ifelse(imp_method %in% c('SVD', 'KNN', 'PPCA', 'BPCA', 'LLS'), 5, NA)),
                        normalization=norm_method, transform=tran_methods)
                ))
                expect_correct_se(
                    se, processed_se, exclude_missing=TRUE,
                    replace_na_method=imp_method, normalization=norm_method, transform=tran_methods)
            }
        }
    }
})


# Test error handling
test_that("data_process handles errors correctly", {
    expect_snapshot(
        data_process(
            se, exclude_missing="TRUE", exclude_missing_pct=70, replace_na_method='min',
            replace_na_method_ref=0.5, normalization='Percentage', transform='log10'),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=101,
            replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10'),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method="invalid_method", replace_na_method_ref=0.5, normalization='Percentage', transform='log10'),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref='invalid_value', normalization='Percentage', transform='log10'),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref=0.5, normalization="invalid_method"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref=0.5, normalization="Percentage", transform="invalid_method"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='BPCA', replace_na_method_ref=11, normalization="Percentage", transform="log10"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='KNN', replace_na_method_ref=11, normalization="Percentage", transform="log10"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='LLS', replace_na_method_ref=27,
            normalization="Percentage", transform="log10"))

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref=1.1, normalization="Percentage", transform="log10"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='PPCA', replace_na_method_ref=11, normalization="Percentage", transform="log10"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='QRILC', replace_na_method_ref=1.1, normalization="Percentage", transform="log10"),
        error=TRUE)

    expect_snapshot(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='SVD', replace_na_method_ref=11, normalization="Percentage", transform="log10"),
        error=TRUE)

})
