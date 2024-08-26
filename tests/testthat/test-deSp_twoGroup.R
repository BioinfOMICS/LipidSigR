library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("de_data_twoGroup")
processed_se <- data_process(de_data_twoGroup, exclude_missing=TRUE,
    exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info = "Object is not an SummarizedExperiment")
    # Check if mandatory components are present
    expect_true("assays" %in% slotNames(se),
                info = "Assays slot is missing")
    # expect_true("rowData" %in% slotNames(se),
    #             info = "rowData slot is missing")
    expect_true("colData" %in% slotNames(se),
                info = "colData slot is missing")
    expect_true("metadata" %in% slotNames(se),
                info = "metadata slot is missing")
}

expect_stat <- function(se){
    # colnames for statistics
    StatCol <- c("feature", "mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp", "FC", "log2FC", "method",
                 "statistic", "pval", "negLog10pval", "padj", "negLog10padj", "sig_pval", "sig_padj")
    # output data frame
    all_deSp_result <- metadata(se)$all_deSp_result
    sig_deSp_result <- metadata(se)$sig_deSp_result
    processed_abundance <- metadata(se)$processed_abundance
    # check se format
    expect_s3_class(all_deSp_result, "data.frame")
    expect_s3_class(sig_deSp_result, "data.frame")
    expect_s3_class(processed_abundance, "data.frame")
    # check expected colnames
    expect_identical(colnames(all_deSp_result), StatCol)
    expect_identical(colnames(sig_deSp_result), StatCol)
}

expect_deSp_plots <- function(se, res){
    all_deSp_result <- metadata(se)$all_deSp_result
    sig_deSp_result <- metadata(se)$sig_deSp_result
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("interactive_de_lipid" %in% names(res), info = "'interactive_de_lipid' should be in the list")
    expect_true("interactive_maPlot" %in% names(res), info = "'interactive_maPlot' should be in the list")
    expect_true("interactive_volcanoPlot" %in% names(res), info = "'interactive_volcanoPlot' should be in the list")
    expect_true("static_de_lipid" %in% names(res), info = "'static_de_lipid' should be in the list")
    expect_true("static_maPlot" %in% names(res), info = "'static_maPlot' should be in the list")
    expect_true("static_volcanoPlot" %in% names(res), info = "'static_volcanoPlot' should be in the list")
    expect_true("table_de_lipid" %in% names(res), info = "'table_de_lipid' should be in the list")
    expect_true("table_ma_volcano" %in% names(res), info = "'table_ma_volcano' should be in the list")

    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_maPlot, "plotly")
    expect_s3_class(res$static_maPlot, "ggplot")
    expect_s3_class(res$table_ma_volcano, "data.frame")
    if (nrow(sig_deSp_result) == 0) {
        expect_true(is.character(res$interactive_de_lipid))
        expect_true(is.character(res$static_de_lipid))
        expect_true(is.character(res$table_de_lipid))
    } else {
        expect_s3_class(res$interactive_de_lipid, "plotly")
        expect_s3_class(res$static_de_lipid, "ggplot")
        expect_s3_class(res$table_de_lipid, "data.frame")
    }
    if (isTRUE("sig_FC" %in% colnames(all_deSp_result))) {
        expect_true(is.character(res$interactive_volcanoPlot))
        expect_true(is.character(res$static_volcanoPlot))
    } else {
        expect_s3_class(res$interactive_volcanoPlot, "plotly")
        expect_s3_class(res$static_volcanoPlot, "ggplot")
    }
}

# Test for basic functionality and output structure
test_that("deSp_twoGroup returns a SummarizedExperiment object with expected metadata", {
    result <- deSp_twoGroup(
        processed_se, ref_group='ctrl', test='t-test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    expect_se_format(result)
    expect_stat(result)
})

# Test for parameter paired
test_that("deSp_twoGroup handles paired samples", {
    sub <- processed_se[, SummarizedExperiment::colData(processed_se)$sample_name %in% c("control_01", "control_02", "hfref_patient_01", "hfref_patient_02")]
    SummarizedExperiment::colData(sub)$pair <- c("1", "1", "2", "2")
    result <- deSp_twoGroup(
        processed_se, ref_group='ctrl', test='t-test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    expect_se_format(result)
    expect_stat(result)
})

# Test for parameter handling transformation types
test_that("deSp_twoGroup handles different normalization and transformation types", {
    transform_type <- c('log10', 'square', 'cube')
    for (i in seq(transform_type)) {
        res <- deSp_twoGroup(
            processed_se, ref_group='ctrl', test='t-test', significant='padj',
            p_cutoff=0.05, FC_cutoff=2, transform=transform_type[i])
        expect_se_format(res)
        expect_stat(res)
        expect_false(identical(assay(res), assay(processed_se) ))
    }
})

# Test for statistical testing and p-value adjustment
test_that("deSp_twoGroup performs statistical testing correctly", {
    result_wilcox <- deSp_twoGroup(
        processed_se, ref_group='ctrl', test='Wilcoxon test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    expect_se_format(result_wilcox)
    expect_stat(result_wilcox)
})

# Test for only 2 samples
test_that("deSp_twoGroup handles only 2 samples", {
    sub <- processed_se[, SummarizedExperiment::colData(processed_se)$sample_name %in% c("control_01", "hfref_patient_01")]
    result <- deSp_twoGroup(
        sub, ref_group='ctrl', test='Wilcoxon test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    expect_se_format(result)
    expect_stat(result)
})

# Test for DE plots
deSp_se <- deSp_twoGroup(
    processed_se, ref_group='ctrl', test='t-test', significant='padj',
    p_cutoff=0.05, FC_cutoff=2, transform='log10')

test_that("plot_deSp_twoGroup returns expected plots", {
    res <- suppressWarnings(
        plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
})

test_that("plot_deSp_twoGroup returns expected plots when significant lipid < 20", {
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][19, ]
    res <- suppressWarnings(
        plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
})

test_that("plot_deSp_twoGroup handles only 1 significant lipid.", {
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, ]
    res <- suppressWarnings(
        plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
})

test_that("plot_deSp_twoGroup handles no significant lipids condition", {
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- data.frame()
    res <- plot_deSp_twoGroup(deSp_se)
    expect_deSp_plots(deSp_se, res)
})
