library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("de_data_twoGroup")
se <- de_data_twoGroup

expect_correct_se <- function(
        input_se, processed_se, exclude_missing, replace_na_method, normalization) {
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
    if (exclude_missing==FALSE && replace_na_method=='none' && normalization=='none') {
        expect_identical(
            SummarizedExperiment::assay(input_se), SummarizedExperiment::assay(processed_se))
    # } else {
    #     expect_false(
    #         identical(SummarizedExperiment::assay(input_se), SummarizedExperiment::assay(processed_se)))
    }
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
})

# Test for all combinations of imputation and normalization methods
test_that("data_process runs successfully on all combinations of imputation and normalization methods", {
    imputation_methods <- c('QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', 'RandomForest')
    #imputation_methods <- c('QRILC', 'SVD', 'KNN', 'min', 'median', 'PPCA', 'BPCA', 'RandomForest')
    #imputation_methods <- 'IRMI'
    normalization_methods <- c('Percentage', 'PQN', 'Quantile', 'Sum', 'Median', 'none')

    for (imp_method in imputation_methods) {
        for (norm_method in normalization_methods) {
            processed_se <- suppressWarnings(suppressMessages(
                data_process(
                    se, exclude_missing=TRUE, exclude_missing_pct=70,
                    replace_na_method=imp_method,
                    replace_na_method_ref=ifelse(
                        imp_method %in% c('QRILC', 'min'), 0.5,
                        ifelse(imp_method %in% c('SVD', 'KNN', 'PPCA', 'BPCA'), 5, NULL)),
                    normalization=norm_method)

            ))
            expect_correct_se(
                se, processed_se, exclude_missing=TRUE,
                replace_na_method=imp_method, normalization=norm_method)
        }
    }

    # for (imp_method in imputation_methods) {
    #     for (norm_method in normalization_methods) {
    #         test_name <- paste("Testing", imp_method, "imputation with", norm_method, "normalization")
    #
    #         test_that(test_name, {
    #             tryCatch({
    #                 processed_se <- suppressWarnings(suppressMessages(
    #                     data_process(
    #                         se, exclude_missing=TRUE, exclude_missing_pct=70,
    #                         replace_na_method=imp_method,
    #                         replace_na_method_ref=ifelse(
    #                             imp_method %in% c('QRILC', 'min'), 0.5,
    #                             ifelse(imp_method %in% c('SVD', 'KNN', 'PPCA', 'BPCA'), 5, NULL)),
    #                         normalization=norm_method)
    #
    #                 ))
    #                 expect_correct_se(
    #                     se, processed_se, exclude_missing=TRUE,
    #                     replace_na_method=imp_method, normalization=norm_method)
    #
    #                 # Additional checks specific to normalization method
    #                 # if (norm_method == 'Percentage') {
    #                 #     expect_true(all(colSums(SummarizedExperiment::assay(processed_se)) %>% near(100)))
    #                 # }
    #             }, error=function(e) {
    #                 fail(paste("Error in combination:", imp_method, "imputation,", norm_method, "normalization:", e$message))
    #             })
    #         })
    #     }
    # }
})

# Test error handling
test_that("data_process handles errors correctly", {
    expect_error(
        data_process(
            se, exclude_missing="TRUE", exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage'),
        "exclude_missing must be a logical value."
    )

    expect_error(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=101,
            replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage'),
        "exclude_missing_pct must be a numeric value between 5 and 100."
    )

    expect_error(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method="invalid_method", replace_na_method_ref=0.5, normalization='Percentage'),
        "replace_na_method must be one of 'none', 'QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', or 'RandomForest'."
    )

    expect_error(
        data_process(
            se, exclude_missing=TRUE, exclude_missing_pct=70,
            replace_na_method='min', replace_na_method_ref=0.5, normalization="invalid_method"),
        "normalization must be one of 'Percentage', 'PQN', 'Quantile', 'Sum', 'Median', or 'none'."
    )
})
