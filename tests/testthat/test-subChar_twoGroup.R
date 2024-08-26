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
# get char
char_list <- c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH')
subChar_list <- list_lipid_char(processed_se)$common_list
names(subChar_list) <- NULL
subChar_list <- subChar_list[!(subChar_list %in% c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH'))]

expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info = "Object is not an SummarizedExperiment")
    expect_true("metadata" %in% slotNames(se),
                info = "metadata slot is missing")
    # Check the meta data
    char <- metadata(se)$char
    subChar <- metadata(se)$subChar
    all_deChar_result <- metadata(se)$all_deChar_result
    sig_deChar_result <- metadata(se)$sig_deChar_result
    processed_abundance <- metadata(se)$processed_abundance
    expect_s3_class(all_deChar_result, "data.frame")
    expect_s3_class(sig_deChar_result, "data.frame")
    expect_s3_class(processed_abundance, "data.frame")
    # col names
    char_col <- c(
        "post_hoc_sig_padj", "sub_characteristic", "sub_feature", "characteristic", "feature", "mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp",
        "FC", "log2FC", "method", "fval_2factors", "pval_2factors", "fval_feature", "pval_feature", "fval_group", "pval_group", "post_hoc_test",
        "post_hoc_statistic", "post_hoc_pval", "post_hoc_negLog10pval", "post_hoc_padj", "post_hoc_negLog10padj", "post_hoc_sig_pval")
    expect_true(all(colnames(all_deChar_result) %in% char_col))
    expect_true(all(colnames(sig_deChar_result) %in% char_col))
}

expect_plots <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("interactive_barPlot" %in% names(res), info = "'interactive_barPlot' should be in the list")
    expect_true("static_barPlot" %in% names(res), info = "'static_barPlot' should be in the list")
    expect_true("interactive_barPlot_sqrt" %in% names(res), info = "'interactive_barPlot_sqrt' should be in the list")
    expect_true("static_barPlot_sqrt" %in% names(res), info = "'static_barPlot_sqrt' should be in the list")
    expect_true("interactive_linePlot" %in% names(res), info = "'interactive_linePlot' should be in the list")
    expect_true("static_linePlot" %in% names(res), info = "'static_linePlot' should be in the list")
    expect_true("interactive_linePlot_sqrt" %in% names(res), info = "'interactive_linePlot_sqrt' should be in the list")
    expect_true("static_linePlot_sqrt" %in% names(res), info = "'static_linePlot_sqrt' should be in the list")
    expect_true("interactive_boxPlot" %in% names(res), info = "'interactive_boxPlot' should be in the list")
    expect_true("table_barPlot" %in% names(res), info = "'table_barPlott' should be in the list")
    expect_true("table_linePlot" %in% names(res), info = "'table_linePlot' should be in the list")
    expect_true("table_boxPlot" %in% names(res), info = "'table_boxPlot' should be in the list")
    expect_true("table_char_index" %in% names(res), info = "'table_char_index' should be in the list")
    expect_true("table_index_stat" %in% names(res), info = "'table_index_stat' should be in the list")
    # Check result type
    expect_s3_class(res$table_barPlot, "data.frame")
    expect_s3_class(res$table_linePlot, "data.frame")
    expect_s3_class(res$table_boxPlot, "data.frame")
    expect_s3_class(res$table_char_index, "data.frame")
    expect_s3_class(res$table_index_stat, "data.frame")

    expect_s3_class(res$interactive_barPlot, "plotly")
    expect_s3_class(res$interactive_barPlot_sqrt, "plotly")
    expect_s3_class(res$interactive_linePlot, "plotly")
    expect_s3_class(res$interactive_linePlot_sqrt, "plotly")
    expect_s3_class(res$interactive_boxPlot, "plotly")
    expect_s3_class(res$static_barPlot, "ggplot")
    expect_s3_class(res$static_barPlot_sqrt, "ggplot")
    expect_s3_class(res$static_linePlot, "ggplot")
    expect_s3_class(res$static_linePlot_sqrt, "ggplot")
    expect_s3_class(res$static_boxPlot, "ggplot")
}

# Test for basic functionality and output structure
test_that("subChar_twoGroup returns correct output", {
    for (i in (seq(1:5))) {
        res <- subChar_twoGroup(
            processed_se, char=char_list[i], subChar=char_list[i+1],
            ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
            FC_cutoff=1, transform='log10')
        expect_se_format(res)
    }
    for (i in seq(char_list) ) {
        for (j in seq(subChar_list)) {
            res <- subChar_twoGroup(
                processed_se, char=char_list[i], subChar=subChar_list[j],
                ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
                FC_cutoff=1, transform='log10')
            expect_se_format(res)
        }
    }
})

test_that("plot_subChar_twoGroup returns correct output", {
    subChar_se <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    subChar_feature_list <- unique(extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)
    for (i in seq(subChar_feature_list)) {
        plot_res <- plot_subChar_twoGroup(subChar_se, subChar_feature=subChar_feature_list[i])
        expect_plots(plot_res)
    }
})
