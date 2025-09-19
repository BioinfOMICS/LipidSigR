library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# load example data
data("de_data_twoGroup")
processed_se <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

data("se_multiGroup")
processed_se_multi <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
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
    expect_true(all(
        c("interactive_barPlot", "static_barPlot", "interactive_barPlot_sqrt", "static_barPlot_sqrt",
          "interactive_linePlot", "static_linePlot", "interactive_linePlot_sqrt", "static_linePlot_sqrt",
          "interactive_boxPlot", "static_boxPlot", "table_barPlot", "table_linePlot", "table_boxPlot",
          "table_char_index", "table_index_stat") %in% names(res)),
        info = "Missing expected elements in the result list")
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
    # Test with different char and subChar combinations
    for (i in seq_along(char_list)) {
        for (j in seq_along(subChar_list)) {
            if (i != j) {
                res <- subChar_twoGroup(
                    processed_se, char=char_list[i], subChar=subChar_list[j],
                    ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
                    FC_cutoff=1, transform='log10')
                expect_se_format(res)
                # test plot
                subChar_feature <- unique(extract_summarized_experiment(res)$all_deChar_result$sub_feature)[1]
                plot_res <- plot_subChar_twoGroup(res, subChar_feature=subChar_feature)
                expect_plots(plot_res)
            }
        }
    }
    # Test with different test, significance, transformations methods
    for(test in c('t-test', 'Wilcoxon test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                res <- subChar_twoGroup(
                    processed_se, char="Total.C", subChar="class",
                    ref_group="ctrl", test=test, significant=significant,
                    p_cutoff=0.05, FC_cutoff=1, transform=transform)
                expect_se_format(res)
            }
        }
    }

    # Test with different p-value cutoffs
    res_p01 <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.01, FC_cutoff=1, transform='log10')
    expect_se_format(res_p01)
    res_p1 <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.1, FC_cutoff=1, transform='log10')
    expect_se_format(res_p1)

    # Test with different FC cutoffs
    res_fc <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1000, transform='log10')
    expect_se_format(res_fc)
})

# Test plot_subChar_twoGroup function
test_that("plot_subChar_twoGroup returns correct output", {
    subChar_se <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    subChar_feature_list <- unique(extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)
    for (feature in subChar_feature_list) {
        plot_res <- plot_subChar_twoGroup(subChar_se, subChar_feature=feature)
        expect_plots(plot_res)
    }
})

# Test error handling
test_that("subChar_twoGroup handles errors correctly", {
    # Test with invalid char
    expect_error(
        subChar_twoGroup(
            processed_se, char="Invalid", subChar="class", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'char' parameter must be a numeric lipid characteristic, such as 'Total.C'. Please execute list_lipid_char function to view the valid lipid characteristics you can input.")

    # Test with invalid subChar
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="Invalid", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'subChar' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input.")

    # Test with same char and subChar
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="Total.C", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'char' and 'subChar' parameters cannot be the same")

    # Test with invalid ref_group
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="class", ref_group="control",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "ref_group must be one of the group names in the 'group' column of the group information table")

    # Test with invalid test method
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="class", ref_group="ctrl",
            test="ttest", significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'test' parameter must be either 't-test', or 'Wilcoxon test'")

    # Test with invalid significant method
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="class", ref_group="ctrl",
            test='t-test', significant="pvalue", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'significant' parameter must be either 'pval' or 'padj'")

    # Test with invalid p_cutoff
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="class", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=2, FC_cutoff=1, transform='log10'),
        "The 'p_cutoff' parameter must be a numeric value between 0 and 1")

    # Test with invalid transform method
    expect_error(
        subChar_twoGroup(
            processed_se, char="Total.C", subChar="class", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform="log2"),
        "The 'transform' parameter must be one of 'none', 'log10', 'square', or 'cube'")
    expect_error(
        subChar_twoGroup(
            processed_se_multi, char="Total.C", subChar="class", ref_group="ctrl",
            test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform="log10"),
        "This function is used for two groups. If your data consists of multi groups, please use deChar_multiGroup functionfor analysis."
    )
})

test_that("plot_subChar_twoGroup handles errors correctly", {
    subChar_se <- subChar_twoGroup(
        processed_se, char="Total.C", subChar="class", ref_group="ctrl",
        test='t-test', significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

    # Test with invalid subChar_feature
    expect_error(
        plot_subChar_twoGroup(subChar_se, subChar_feature="Invalid"),
        "The 'subChar_feature' parameter must be one of the valid subChar features.")
})
