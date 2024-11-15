library(testthat)
library(SummarizedExperiment)
library(tidyverse)
# Load data and create example data sets
data("ml_data")
processed_se <- data_process(
    ml_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info="Object is not an SummarizedExperiment")
    # Check if mandatory components are present
    expect_true("assays" %in% slotNames(se), info="Assays slot is missing")
    # expect_true("rowData" %in% slotNames(se),
    #             info="rowData slot is missing")
    expect_true("colData" %in% slotNames(se), info="colData slot is missing")
    expect_true("metadata" %in% slotNames(se), info="metadata slot is missing")
    ## check metadata
    expect_equal(length(se@metadata), 15)
    var_list <- c('char', "transform", "ranking_method", "ml_method")
    result_list <- c(
        "model_result", "confusion_matrix", "roc_result", "pr_result",
        "mean_roc_result", "feature_importance_result")
    list_res <- c("selected_features", "best_model", "best_model_feature")
    expect_identical(
        names(se@metadata),
        c(var_list, "nfold", "feature_option", result_list, list_res))
    all_valid <- all(
        vapply(var_list, function(res) is.character(S4Vectors::metadata(se)[[res]]), logical(1)),
        vapply(result_list, function(res) is.data.frame(S4Vectors::metadata(se)[[res]]), logical(1)),
        vapply(list_res, function(res) is.list(S4Vectors::metadata(se)[[res]]), logical(1))
    )
    expect_true(all_valid)
    expect_equal(length(S4Vectors::metadata(se)[["feature_option"]]), 7)
}

test_that("Input parameter validation works correctly", {
    # Test invalid char parameter
    expect_error(
        ml_model(
            processed_se, char='char', transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    expect_error(
        ml_model(
            processed_se, char=NULL, transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    # Test invalid transform parameter
    expect_error(
        ml_model(
            processed_se, char="none", transform=NULL, ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL),
        "transform must be one of 'none', 'log10', 'square', or 'cube'."
    )
    # Test invalid split_prop parameter
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0, nfold=3, alpha=NULL),
        "split_prop must be a numeric value between 0.1 and 0.5."
    )
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=1, nfold=3, alpha=NULL),
        "split_prop must be a numeric value between 0.1 and 0.5."
    )
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop='0.3', nfold=3, alpha=NULL),
        "split_prop must be a numeric value between 0.1 and 0.5."
    )
    # Test invalid nfold parameter
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=0, alpha=NULL),
        "nfold must be a numeric value >= 1."
    )
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold='3', alpha=NULL),
        "nfold must be a numeric value >= 1."
    )
    # Test invalid alpha parameter with ElasticNet
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=NULL),
        "alpha must be a numeric value between 0 and 1."
    )
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=-1),
        "alpha must be a numeric value between 0 and 1."
    )
    expect_error(
        ml_model(
            processed_se, char="none", transform='log10', ranking_method='Random_forest',
            ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=2),
        "alpha must be a numeric value between 0 and 1."
    )
})

test_that("Ranking methods & machine learning methods work correctly for all conditions", {
    # Test all ML methods
    ml_methods <- c('Random_forest', 'SVM', 'Lasso', 'Ridge', 'xgboost', 'ElasticNet')
    ranking_methods <- c('p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet')
    for (rank in ranking_methods) {
        for(method in ml_methods) {
            res <- suppressWarnings(
                ml_model(
                    processed_se,
                    char=c("class", "Total.OH", "Total.DB", "FA.OH", "FA.DB", "FA.C"),
                    ranking_method=rank,
                    ml_method=method, split_prop=0.3, nfold=3, alpha=0.5)
            )
            expect_se_format(res)
            expect_equal(S4Vectors::metadata(res)$ranking_method, rank)
            expect_equal(S4Vectors::metadata(res)$ml_method, method)
        }
    }
    # Test invalid ranking method
    expect_error(
        ml_model(
            processed_se, char='none', transform='log10', ranking_method=NULL,
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL),
        "ranking_method must be one of 'p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet'."
    )
    expect_error(
        ml_model(
            processed_se, char='none', transform='log10', ranking_method='ElasticNet',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha='1'),
        "alpha must be a numeric value between 0 and 1."
    )
    # Test invalid ML method
    expect_error(
        ml_model(
            processed_se, char='none', transform='log10', ranking_method='Random_forest',
            ml_method=NULL, split_prop=0.3, nfold=3, alpha=NULL),
        "ml_method must be one of 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost'."
    )
})

test_that("Data transformation options work correctly", {
    # Test all transform options
    transform_methods <- c('none', 'log10', 'square', 'cube')
    for(method in transform_methods) {
        res <- ml_model(
            processed_se, char='none', transform=method, ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL)
        expect_se_format(res)
        expect_equal(S4Vectors::metadata(res)$transform, method)
    }
    # Test invalid transform method
    expect_error(
        ml_model(
            processed_se, char='none', transform=NULL, ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL),
        "transform must be one of 'none', 'log10', 'square', or 'cube'."
    )
})

test_that("Cross-validation parameters work correctly", {
    # Test different split proportions
    split_props <- c(0.1, 0.3, 0.5)
    for(prop in split_props) {
        res <- ml_model(
            processed_se, char='none', transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=prop, nfold=3)
        expect_se_format(res)
    }
    # Test different number of folds
    nfolds <- c(3, 5, 10)
    for(fold in nfolds) {
        res <- ml_model(
            processed_se, char='none', transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=fold)
        expect_se_format(res)
        expect_equal(length(S4Vectors::metadata(res)$selected_features), fold)
    }
})

test_that("ElasticNet alpha parameter works correctly", {
    # Test valid alpha values
    alpha_values <- c(0.2, 0.5, 0.8)
    for(alpha in alpha_values) {
        res <- suppressWarnings(
            ml_model(
                processed_se, char='none', transform='log10', ranking_method='Random_forest',
                ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=alpha)
        )
        expect_se_format(res)
    }
    # Test boundary values
    res_lasso <- suppressWarnings(
        ml_model(
            processed_se, char='none', transform='log10', ranking_method='Random_forest',
            ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=1)
    )
    res_ridge <- suppressWarnings(
        ml_model(
            processed_se, char='none', transform='log10', ranking_method='Random_forest',
            ml_method='ElasticNet', split_prop=0.3, nfold=3, alpha=0)
    )
    expect_se_format(res_lasso)
    expect_se_format(res_ridge)
})

test_that("All char can works correctly.", {
    char_list <- c("class", "Total.OH", "Total.DB", "FA.OH", "FA.DB", "FA.C")
    for (char in char_list) {
        res <- suppressWarnings(
            ml_model(
                processed_se, char=char, transform='log10', ranking_method='Random_forest',
                ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL)
        )
        expect_se_format(res)
    }
    res <- suppressWarnings(
        ml_model(
            processed_se, char=char_list, transform='log10', ranking_method='Random_forest',
            ml_method='Random_forest', split_prop=0.3, nfold=3, alpha=NULL)
    )
    expect_se_format(res)
})
