library(testthat)
library(SummarizedExperiment)
library(tidyverse)
# Load data and create example data sets
data(ml_data)
processed_se <- data_process(
    ml_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
ml_se <- ml_model(
    processed_se, char=c("class","Total.DB"), transform='log10',
    ranking_method='Random_forest', ml_method='Random_forest', split_prop=0.3,
    nfold=3, alpha=NULL)

expect_pr_plots <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_mean_auc", "static_mean_auc", "initeractive_pr",
          "static_pr", "table_mean_auc_plot", "table_pr"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_mean_auc, "plotly")
    expect_s3_class(res$static_mean_auc, "ggplot")
    expect_s3_class(res$table_mean_auc_plot, "data.frame")
    expect_s3_class(res$initeractive_pr, "plotly")
    expect_s3_class(res$static_pr, "ggplot")
    expect_s3_class(res$table_pr, "data.frame")
}

test_that("plot_ml_pr returns expected plots", {
    feature_opt <- S4Vectors::metadata(ml_se)$feature_option
    for (feature_num in feature_opt) {
        res <- plot_ml_pr(ml_se, feature_num)
        expect_pr_plots(res)
    }
})

test_that("plot_ml_pr can handle wrong input", {
    expect_error(
        plot_ml_pr(ml_se, feature_num=NULL),
        "feature_num must be a numeric value chosen from the following options: 2, 3, 5, 10, 20, 50, 100."
    )
})
