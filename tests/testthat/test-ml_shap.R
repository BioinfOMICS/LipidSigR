library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(SHAPforxgboost)
# Load data and create example data sets
data("ml_data")
processed_se <- data_process(
    ml_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
ml_se <- ml_model(
    processed_se, char='none', transform='log10',
    ranking_method='Random_forest', ml_method='Random_forest', split_prop=0.3,
    nfold=3, alpha=NULL)

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
    expect_equal(length(se@metadata), 4)
    char_list <- c("feature_num", "nsim")
    df_list <- c("shap_result", "shap_score")
    expect_identical(names(se@metadata), c(char_list, df_list))
    all_valid <- all(
        vapply(char_list, function(res) is.numeric(S4Vectors::metadata(se)[[res]]), logical(1)),
        vapply(df_list, function(res) is.data.frame(S4Vectors::metadata(se)[[res]]), logical(1))
    )
    expect_true(all_valid)
}

expect_shap_plots <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_feature_importance", "static_feature_importance", "interactive_summary_plot",
          "static_summary_plot", "table_feature_importance", "table_summary_plot"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_feature_importance, "plotly")
    expect_s3_class(res$static_feature_importance, "ggplot")
    expect_s3_class(res$interactive_summary_plot, "plotly")
    expect_s3_class(res$static_summary_plot, "ggplot")
    expect_s3_class(res$table_feature_importance, "data.frame")
    expect_s3_class(res$table_summary_plot, "data.frame")
}

expect_shap_sample <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_sample_feature_importance", "static_sample_feature_importance",
          "table_sample_feature_importance"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_sample_feature_importance, "plotly")
    expect_s3_class(res$static_sample_feature_importance, "ggplot")
    expect_s3_class(res$table_sample_feature_importance, "data.frame")
}

expect_shap_force <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_forcePlot", "static_forcePlot", "table_forcePlot"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_forcePlot, "plotly")
    expect_s3_class(res$static_forcePlot, "ggplot")
    expect_s3_class(res$table_forcePlot, "data.frame")
}

expect_shap_dependence <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_dependence_plot", "static_dependence_plot", "table_dependence_plot"))
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_dependence_plot, "plotly")
    expect_s3_class(res$static_dependence_plot, "ggplot")
    expect_s3_class(res$table_dependence_plot, "data.frame")
}

# Start tests
test_that("ml_shap & plot_ml_shap works correctly and can handle error inputs", {
    ml_methods <- c('Random_forest', 'SVM', 'Lasso', 'Ridge', 'xgboost', 'ElasticNet')
    for(method in ml_methods) {
        ml_res <- suppressWarnings(
            ml_model(
                processed_se,
                char=c("class", "Total.OH", "Total.DB", "FA.OH", "FA.DB", "FA.C"),
                ranking_method='Random_forest',
                ml_method=method, split_prop=0.3, nfold=3, alpha=0.5)
        )
        feature_num <- S4Vectors::metadata(ml_res)$feature_option
        for (num in feature_num) {
            shap_se <- ml_shap(ml_res, feature_num=num, nsim=5)
            expect_se_format(shap_se)
            res <- suppressWarnings( plot_ml_shap(shap_se) )
            expect_shap_plots(res)
        }
    }
    expect_error(
        ml_shap(ml_se, feature_num='200', nsim=5),
        "feature_num must be a numeric value chosen from the following options: "
    )
    expect_error(
        ml_shap(ml_se, feature_num=10, nsim=-1),
        "nsim must be a positive integer."
    )
})

shap_se <- ml_shap(ml_se, feature_num=10, nsim=5)

test_that("plot_shap_sample can return correct plots.", {
    res <- plot_shap_sample(shap_se, sample_id=10)
    expect_shap_sample(res)
    expect_error(
        plot_shap_sample(shap_se, sample_id='10'),
        "sample_id must be one of the ID number in 'ID' column of shap_result."
    )
})

test_that("plot_shap_force can return correct plots.", {
    cluster_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
    for (cluster in cluster_methods) {
        res <- plot_shap_force(shap_se, top_feature=10, cluster_method=cluster, group_num=10)
        expect_shap_force(res)
    }
    expect_error(
        plot_shap_force(shap_se, top_feature=10, cluster_method='cluster', group_num=10),
        'cluster_method must be one of the strings "ward.D", "ward.D2",
         "single", "complete", "average", "mcquitty", "median", "centroid".'
    )
    expect_error(
        plot_shap_force(shap_se, top_feature=10, cluster_method='ward.D', group_num='10'),
        "group_num must be a numeric value between 0 and "
    )
    expect_error(
        plot_shap_force(shap_se, top_feature='10', cluster_method='ward.D', group_num=10),
        "top_feature must be a positive value."
    )
    expect_warning(
        plot_shap_force(shap_se, top_feature=20, cluster_method='ward.D', group_num=10),
        "top_feature is set to greater than 10, only the top 10 samples will be displayed"
    )

})

test_that("plot_shap_dependence can return correct plots.", {
    selected_feature <- as.character(unique(S4Vectors::metadata(shap_se)$shap_result$variable))
    res <- suppressWarnings(
        plot_shap_dependence(
            shap_se, feature=selected_feature[1], shap_feature=selected_feature[2],
            interaction_index=selected_feature[2])
    )
    expect_shap_dependence(res)
    expect_error(
        plot_shap_dependence(
            shap_se, feature='feature1', shap_feature=selected_feature[2],
            interaction_index=selected_feature[2]),
        "feature must be one of the feature names in the 'variable' column of shap_result."
    )
    expect_error(
        plot_shap_dependence(
            shap_se, feature=selected_feature[1], shap_feature=NULL,
            interaction_index=selected_feature[2]),
        "shap_feature must be one of the feature names in the 'variable' column of shap_result."
    )
    expect_error(
        plot_shap_dependence(
            shap_se, feature=selected_feature[1], shap_feature=selected_feature[2],
            interaction_index='abd'),
        "interaction_index must be one of the feature names in the 'variable' column of shap_result."
    )
})
