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

expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info = "Object is not an SummarizedExperiment")
    expect_true("metadata" %in% slotNames(se),
                info = "metadata slot is missing")
    # Check the meta data
    char <- metadata(se)$char
    all_deChar_result <- metadata(se)$all_deChar_result
    sig_deChar_result <- metadata(se)$sig_deChar_result
    processed_abundance <- metadata(se)$processed_abund
    expect_s3_class(all_deChar_result, "data.frame")
    expect_s3_class(sig_deChar_result, "data.frame")
    expect_s3_class(processed_abundance, "data.frame")
    # col names
    # char_col <- c(
    #     "feature", "mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp", "FC", "log2FC",
    #     "method", "anova.fval", "anova.pval", "post.hoc.method", "post.hoc.pval",
    #     "post.hoc.mlog10p", "post.hoc.padj", "post.hoc.mlog10padj", "post.hoc.pval.signif", "post.hoc.padj.signif")
    # bar_col <- c("Category", "Significant", "Mean", "SD", "Group", "max_error_bar", "post_hoc_pvalue", "pvalue_text")
    # expect_identical(colnames(all_deChar_result), char_col)
    # expect_identical(colnames(barTab_sig), bar_col)
}

test_that("deChar_twoGroup function outputs correct structure", {
    for (i in seq(char_list) ) {
        res <- deChar_twoGroup(
            processed_se, char=char_list[i], ref_group="ctrl", test='t-test',
            significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='none')
        expect_se_format(res)
    }
})

# Test for parameter paired
test_that("deChar_twoGroup handles paired samples", {
    sub <- processed_se[, SummarizedExperiment::colData(processed_se)$sample_name %in%
                  c("control_01", "control_02", "control_03", "control_04", "control_05", "control_06", "control_07", "control_08", "control_09", "control_10",
                    "hfref_patient_01", "hfref_patient_02", "hfref_patient_03", "hfref_patient_04", "hfref_patient_05", "hfref_patient_06", "hfref_patient_07", "hfref_patient_08", "hfref_patient_09", "hfref_patient_10", "hfref_patient_11", "hfref_patient_12")]
    SummarizedExperiment::colData(sub)$pair <- c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6", "7", "7", "8", "8", "9", "9", "10", "10", "11", "11")
    result <- deChar_twoGroup(
        sub, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    expect_se_format(result)
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
    result_wilcox <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='Wilcoxon test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    expect_se_format(result_wilcox)
    result_p <- deChar_twoGroup(
        processed_se, char="Total.C", ref_group="ctrl", test='Wilcoxon test',
        significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    expect_se_format(result_p)
})

# Test for DE_char plotting
# test_that("deChar_plot performs output plots correctly", {
#     de_char_se <- deChar(char_data, char="Total.C")
#     char_plot <- deChar_plot(de_char_se, print_static=TRUE)
#     de_char_class <- deChar(char_data, char="class")
#     char_plot_S <- deChar_plot(de_char_class)
#     de_char_p <- deChar(char_data, sig_stat='p', char="Total.C")
#     char_plot <- deChar_plot(de_char_p, print_static=TRUE)
# })


