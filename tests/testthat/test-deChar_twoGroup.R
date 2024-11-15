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
char_list <- list_lipid_char(processed_se)$deChar_list
names(char_list) <- NULL

data("se_multiGroup")
processed_se_multi <- data_process(
    se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')


expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info = "Object is not an SummarizedExperiment")
    expect_true("metadata" %in% slotNames(se),
                info = "metadata slot is missing")
    # Check metadata contents
    expect_type(metadata(se)$char, "character")
    expect_s3_class(metadata(se)$all_deChar_result, "data.frame")
    expect_true(is.data.frame(metadata(se)$sig_deChar_result) || is.character(metadata(se)$sig_deChar_result))
    expect_s3_class(metadata(se)$processed_abundance, "data.frame")

    # Check for required columns in all_deChar_result
    required_cols <- c("feature", "mean_ctrl", "mean_exp", "FC", "log2FC", "post_hoc_pval", "post_hoc_padj")
    expect_true(all(required_cols %in% colnames(metadata(se)$all_deChar_result)))
}

expect_plot_result <- function(plot_result, char){
    expect_type(plot_result, "list")
    if (char %in% c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH')) {
        expect_named(plot_result,
                     c("static_barPlot", "static_barPlot_sqrt", "static_linePlot",
                       "static_linePlot_sqrt", "static_boxPlot", "interactive_barPlot",
                       "interactive_barPlot_sqrt", "interactive_linePlot",
                       "interactive_linePlot_sqrt", "interactive_boxPlot",
                       "table_barPlot", "table_linePlot", "table_boxPlot",
                       "table_char_index", "table_index_stat"))
        expect_s3_class(plot_result$static_linePlot, "ggplot")
        expect_s3_class(plot_result$static_linePlot_sqrt, "ggplot")
        expect_s3_class(plot_result$static_boxPlot, "ggplot")

        expect_s3_class(plot_result$interactive_linePlot, "plotly")
        expect_s3_class(plot_result$interactive_linePlot_sqrt, "plotly")
        expect_s3_class(plot_result$interactive_boxPlot, "plotly")

        expect_s3_class(plot_result$table_linePlot, "data.frame")
        expect_s3_class(plot_result$table_boxPlot, "data.frame")

    } else {
        expect_named(plot_result,
                     c("static_barPlot", "static_barPlot_sqrt",
                       "interactive_barPlot", "interactive_barPlot_sqrt",
                       "table_barPlot"))
    }
    expect_s3_class(plot_result$static_barPlot, "ggplot")
    expect_s3_class(plot_result$static_barPlot_sqrt, "ggplot")
    expect_s3_class(plot_result$interactive_barPlot, "plotly")
    expect_s3_class(plot_result$interactive_barPlot_sqrt, "plotly")
    expect_s3_class(plot_result$table_barPlot, "data.frame")
}


test_that("deChar_twoGroup function outputs correct structure for all characteristics", {
    for (char in char_list) {
        res <- suppressWarnings(
            deChar_twoGroup(
                processed_se, char=char, ref_group="ctrl", test='t-test',
                significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='none')
        )
        expect_se_format(res)
        expect_identical(metadata(res)$char, char)
    }
})

# Test for parameter paired
test_that("deChar_twoGroup handles paired samples correctly", {
    sub <- processed_se[, colData(processed_se)$sample_name %in%
                            c("control_01", "control_02", "control_03", "control_04", "control_05",
                              "hfref_patient_01", "hfref_patient_02", "hfref_patient_03", "hfref_patient_04", "hfref_patient_05")]
    colData(sub)$pair <- rep(1:5, 2)
    result <- suppressWarnings(
        deChar_twoGroup(
            sub, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    )
    expect_se_format(result)
    expect_true("pair" %in% colnames(colData(result)))
})

# Test for parameter handling (transformation types)
test_that("deChar_twoGroup handles different transformation types", {
    transform_type <- c('none', 'log10', 'square', 'cube')
    for (i in seq(transform_type)) {
        res <- deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform=transform_type[i])
        expect_se_format(res)
    }
})

# Test for statistical testing and p-value adjustment
test_that("deChar_twoGroup performs statistical testing correctly", {
    result_ttest <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    result_wilcox <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='Wilcoxon test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')

    expect_se_format(result_ttest)
    expect_se_format(result_wilcox)
    expect_false(identical(metadata(result_ttest)$all_deChar_result$post_hoc_pval,
                           metadata(result_wilcox)$all_deChar_result$post_hoc_pval))
})

test_that("deChar_twoGroup handles different significance methods", {
    result_padj <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    result_pval <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

    expect_se_format(result_padj)
    expect_se_format(result_pval)
    expect_false(identical(metadata(result_padj)$sig_deChar_result,
                           metadata(result_pval)$sig_deChar_result))
})

test_that("deChar_twoGroup handles edge cases", {
    # Test with no significant results
    res_no_sig <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=0.000001, FC_cutoff=100, transform='log10')
    expect_se_format(res_no_sig)
    expect_equal(metadata(res_no_sig)$sig_deChar_result, "No significant lipids")

    # Test with all significant results
    res_all_sig <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=1, FC_cutoff=1, transform='log10')
    expect_se_format(res_all_sig)
    expect_equal(nrow(metadata(res_all_sig)$sig_deChar_result),
                 nrow(metadata(res_all_sig)$all_deChar_result))
})

test_that("deChar_twoGroup handles invalid inputs", {
    expect_error(
        deChar_twoGroup(
            processed_se, char="Invalid", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "The 'char' parameter is not in the list of lipid characteristics")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="invalid", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "ref_group must be one of the group names")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test="invalid",
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "test must be one of 't-test' or 'Wilcoxon test'")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test='t-test',
            significant="invalid", p_cutoff=0.05, FC_cutoff=1, transform='log10'),
                 "significant must be one of 'pval' or 'padj'")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=2, FC_cutoff=1, transform='log10'),
                 "p_cutoff must be a numeric value between 0 and 1")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=0),
        "FC_cutoff must be a numeric value >= 1.")
    expect_error(
        deChar_twoGroup(
            processed_se, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=1, FC_cutoff=1, transform='log2'),
        "The 'char' is not a ratio-based characteristic. The 'transform' parameter must be one of 'none', 'log10', 'square', or 'cube'."
    )
    expect_error(
        deChar_twoGroup(
            processed_se, char="Chains odd/even ratio", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=1, FC_cutoff=1, transform='log10'),
        "The 'char' is a ratio-based characteristic. The 'transform' parameter must be either 'none' or 'log2'."
    )
    expect_error(
        deChar_twoGroup(
            processed_se_multi, char="Total.C", ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=1, FC_cutoff=1, transform='log10'),
        "This function is used for two groups. If your data consists of multi groups, please use deChar_nultiGroup function for analysis."
    )
})

## test plot_deChar_twoGroup
test_that("plot_deChar_twoGroup function produces correct output", {
    deChar_result <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=1, FC_cutoff=1, transform='log10')
    plot_result <- plot_deChar_twoGroup(deChar_result)
    expect_plot_result(plot_result, char="Total.C")
})

test_that("plot_deChar_twoGroup handles non-numeric characteristics", {
    deChar_result <- deChar_twoGroup(
        processed_se, char="class", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=1, FC_cutoff=1, transform='log10')
    expect_message(
        plot_result <- plot_deChar_twoGroup(deChar_result),
        "The 'char' parameter is not a numeric characteristic so it is unable to perform line plot and box plot."
    )
    expect_plot_result(plot_result, char="class")
})

test_that("plot_deChar_twoGroup handles invalid inputs", {
    deChar_se <- deChar_multiGroup(
        processed_se_multi, char='class', ref_group='ctrl', post_hoc='One-way ANOVA',
        post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10')
    expect_error(
        plot_deChar_twoGroup(deChar_se),
        "This function is used for two groups. If your data consists of two groups, please use plot_deChar_multiGroup function for analysis."
    )
})
