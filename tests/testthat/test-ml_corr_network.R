library(testthat)
library(SummarizedExperiment)
library(tidyverse)

# Load data and create example data sets
data("ml_data")
processed_se <- data_process(
    ml_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
ml_se <- ml_model(
    processed_se, char='none', transform='log10',
    ranking_method='Random_forest', ml_method='Random_forest', split_prop=0.3,
    nfold=3, alpha=NULL)

expect_network <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_correlation_network", "static_correlation_network",
          "edge_table", "node_table"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_correlation_network, "visNetwork")
    expect_s3_class(res$static_correlation_network, "ggplot")
    expect_s3_class(res$edge_table, "data.frame")
    expect_s3_class(res$node_table, "data.frame")
}

# Start tests
test_that("ml_corr_network functions works", {
    res <- ml_corr_network(
        ml_se, feature_importance='Algorithm-based', correlation='pearson',
        edge_cutoff=0, feature_num=10, nsim=5)
    expect_network(res)
})

test_that("ml_corr_network can handle input errors.", {
    expect_error(
        ml_corr_network(
            ml_se, feature_importance='correlation', correlation='pearson',
            edge_cutoff=0, feature_num=10, nsim=5),
        'feature_importance must be "Algorithm-based" or "SHAP".'
    )
    expect_error(
        ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation=NULL,
            edge_cutoff=0, feature_num=10, nsim=5),
        'correlation must be one of "pearson", "kendall", or "spearman".'
    )
    expect_error(
        ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation='pearson',
            edge_cutoff=0, feature_num='10', nsim=5),
        'feature_num must be a numeric value chosen from the following options: 2, 3, 5, 10, 20, 50, 88.'
    )
    expect_error(
        ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation='pearson',
            edge_cutoff=0, feature_num=10, nsim=-1),
        'nsim must be a integer >= 0.'
    )
    expect_error(
        ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation='pearson',
            edge_cutoff='0', feature_num=10, nsim=5),
        'edge_cutoff must be a numeric value between 0 and 1.'
    )
})

test_that("ml_corr_network can handle all feature_importance methods and correlation methods.", {
    res <- ml_corr_network(
        ml_se, feature_importance='SHAP', correlation='pearson',
        edge_cutoff=0, feature_num=10, nsim=5)
    expect_network(res)
    ml_method_list <- c('SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost')
    for (method in ml_method_list) {
        ml_se <- suppressWarnings(
            ml_model(
                processed_se, char=c("class","Total.DB"), transform='log10',
                ranking_method='Random_forest', ml_method=method,
                split_prop=0.3, nfold=5, alpha=0.5)
        )
        res <- ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation='pearson',
            edge_cutoff=0, feature_num=10, nsim=5)
        expect_network(res)
    }
    corr_methods <- c('pearson', 'kendall', 'spearman')
    for (method in corr_methods) {
        res <- ml_corr_network(
            ml_se, feature_importance='Algorithm-based', correlation=method,
            edge_cutoff=0, feature_num=10, nsim=5)
        expect_network(res)
    }
})
