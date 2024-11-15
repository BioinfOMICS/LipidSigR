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

expect_eva_plots <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_evaluation_plot", "static_evaluation_plot",
          "table_evaluation", "table_evaluation_plot"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_evaluation_plot, "plotly")
    expect_s3_class(res$static_evaluation_plot, "ggplot")
    expect_s3_class(res$table_evaluation, "data.frame")
    expect_s3_class(res$table_evaluation_plot, "data.frame")
}

test_that("plot_ml_evaluation returns expected plots", {
    eval_method_option <- c(
        'Accuracy', 'Sensitivity', 'Specificity', 'Pos Pred Value',
        'Neg Pred Value', 'Precision', 'Recall', 'F1', 'Prevalence',
        'Detection Rate', 'Detection Prevalence', 'Balanced Accuracy')
    for (eval_method in eval_method_option) {
        res <- plot_ml_evaluation(ml_se, eval_method)
        expect_eva_plots(res)
    }
})


test_that("plot_ml_evaluation handles wrong input", {
    expect_error(
        plot_ml_evaluation(ml_se, eval_method=NULL),
        'eval_method must be one of the strings "Accuracy", "Sensitivity",
         "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision",
         "Recall", "F1", "Prevalence", "Detection Rate",
         "Detection Prevalence", or "Balanced Accuracy".'
    )
})



