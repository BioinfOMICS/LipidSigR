library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(tidyr)
library(magrittr)

data("de_data_twoGroup")
processed_se <- data_process(de_data_twoGroup, exclude_missing=TRUE,
    exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

data("se_multiGroup")
processed_se_multi <- data_process(
    se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

expect_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info="Object is not an SummarizedExperiment")
    # Check if mandatory components are present
    expect_true("assays" %in% slotNames(se),
                info="Assays slot is missing")
    # expect_true("rowData" %in% slotNames(se),
    #             info="rowData slot is missing")
    expect_true("colData" %in% slotNames(se),
                info="colData slot is missing")
    expect_true("metadata" %in% slotNames(se),
                info="metadata slot is missing")
}

expect_stat <- function(se){
    # colnames for statistics
    StatCol <- c("feature", "mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp", "FC", "log2FC", "method",
                 "statistic", "pval", "negLog10pval", "padj", "negLog10padj", "sig_pval", "sig_padj", "sig_FC")
    # output data frame
    all_deSp_result <- S4Vectors::metadata(se)$all_deSp_result
    sig_deSp_result <- S4Vectors::metadata(se)$sig_deSp_result
    processed_abundance <- S4Vectors::metadata(se)$processed_abundance
    # check se format
    expect_s3_class(all_deSp_result, "data.frame")
    expect_s3_class(processed_abundance, "data.frame")
    # check expected colnames
    expect_true(all(colnames(all_deSp_result) %in% StatCol))
    if (!is.character(sig_deSp_result)) {
        expect_true(all(colnames(sig_deSp_result) %in% StatCol))
        expect_s3_class(sig_deSp_result, "data.frame")
    }
}

expect_deSp_plots <- function(se, res){
    all_deSp_result <- S4Vectors::metadata(se)$all_deSp_result
    sig_deSp_result <- S4Vectors::metadata(se)$sig_deSp_result
    # Check if `res` is a list
    expect_true(is.list(res), info="Output should be a list")
    expect_named(
        res,
        c("interactive_de_lipid", "interactive_maPlot", "interactive_volcanoPlot",
          "static_de_lipid", "static_maPlot", "static_volcanoPlot",
          "table_de_lipid", "table_ma_volcano"))

    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_maPlot, "plotly")
    expect_s3_class(res$static_maPlot, "ggplot")
    expect_s3_class(res$table_ma_volcano, "data.frame")
    if (is.character(sig_deSp_result)) {
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
    expect_equal(S4Vectors::metadata(result)$significant, "padj")
    expect_equal(S4Vectors::metadata(result)$p_cutoff, 0.05)
    expect_equal(S4Vectors::metadata(result)$FC_cutoff, 2)
    expect_equal(S4Vectors::metadata(result)$transform, "log10")
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

test_that("deSp_twoGroup handles different parameter combinations", {
    param_combinations <- tidyr::crossing(
        test = c('t-test', 'Wilcoxon test'),
        significant = c('pval', 'padj'),
        transform = c('none', 'log10', 'square', 'cube')
    )

    for (i in 1:nrow(param_combinations)) {
        result <- deSp_twoGroup(
            processed_se, ref_group='ctrl',
            test=param_combinations$test[i],
            significant=param_combinations$significant[i],
            p_cutoff=0.05, FC_cutoff=2,
            transform=param_combinations$transform[i])
        expect_se_format(result)
        expect_stat(result)
    }
})

test_that("deSp_twoGroup handles edge cases correctly", {
    # Test with only two samples
    sub_se <- processed_se[, SummarizedExperiment::colData(processed_se)$sample_name %in% c("control_01", "hfref_patient_01")]
    result <- deSp_twoGroup(
        sub_se, ref_group='ctrl', test='t-test',
        significant='padj', p_cutoff=0.05, FC_cutoff=1, transform='log10')
    expect_se_format(result)
    expect_stat(result)
    if (!is.character(S4Vectors::metadata(result)$sig_deSp_result)) {
        expect_true("sig_FC" %in% colnames(S4Vectors::metadata(result)$sig_deSp_result))
    }

    # Test with very high FC_cutoff
    result_high_fc <- suppressWarnings(
        deSp_twoGroup(
            processed_se, ref_group='ctrl', test='t-test',
            significant='padj', p_cutoff=0.05, FC_cutoff=1000,
            transform='log10')
    )
    expect_warning(
        deSp_twoGroup(
            processed_se, ref_group='ctrl', test='t-test',
            significant='padj', p_cutoff=0.05, FC_cutoff=1000,
            transform='log10'),
        "This case does not include any significant lipids, which will prevent some subsequent analyses from being performed."
    )
    expect_se_format(result_high_fc)
    expect_stat(result_high_fc)
})

test_that("deSp_twoGroup handles errors correctly", {
    expect_error(deSp_twoGroup(processed_se, ref_group='nonexistent'),
                 "ref_group must be one of the group names in the 'group' column of the group information table.")
    expect_error(deSp_twoGroup(processed_se, ref_group='ctrl', test='invalid_test'),
                 "test must be one of 't-test' or 'Wilcoxon test'.")
    expect_error(deSp_twoGroup(processed_se, ref_group='ctrl', significant='invalid'),
                 "significant must be one of 'pval' or 'padj'.")
    expect_error(deSp_twoGroup(processed_se, ref_group='ctrl', p_cutoff=2),
                 "p_cutoff must be a numeric value between 0 and 1.")
    expect_error(deSp_twoGroup(processed_se, ref_group='ctrl', FC_cutoff=0),
                 "FC_cutoff must be a numeric value >= 1.")
    expect_error(deSp_twoGroup(processed_se, ref_group='ctrl', transform='invalid'),
                 "transform must be one of 'none', 'log10', 'square', or 'cube'.")
    expect_error(
        deSp_twoGroup(
            processed_se_multi, ref_group='ctrl', test='t-test', significant='padj',
            p_cutoff=0.05, FC_cutoff=1, transform='log10'),
        "This function is used for two groups. If your data consists of multi groups, please use deSp_multiGroup function for analysis."

    )
})


# Test for DE plots
deSp_se <- deSp_twoGroup(
    processed_se, ref_group='ctrl', test='t-test', significant='padj',
    p_cutoff=0.05, FC_cutoff=1, transform='log10')

test_that("plot_deSp_twoGroup returns expected plots", {
    res <- suppressWarnings(
        plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
})

test_that("plot_deSp_twoGroup handles different numbers of significant lipids", {
    deSp_se <- deSp_twoGroup(
        processed_se, ref_group='ctrl', test='t-test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')

    # Test with > 20 significant lipids
    res <- suppressWarnings(plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
    expect_lte(nrow(res$table_de_lipid), 20)

    # Test with < 20 significant lipids
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1:19, ]
    res <- suppressWarnings(plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
    expect_equal(nrow(res$table_de_lipid), 19)

    # Test with 1 significant lipid
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, ]
    res <- suppressWarnings(plot_deSp_twoGroup(deSp_se))
    expect_deSp_plots(deSp_se, res)
    expect_equal(nrow(res$table_de_lipid), 1)

    # Test with no significant lipids
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- "No significant lipids."
    res <- plot_deSp_twoGroup(deSp_se)
    expect_deSp_plots(deSp_se, res)
    expect_true(is.character(res$table_de_lipid))
})

test_that("plot_deSp_twoGroup handles edge cases correctly", {
    # Test with only two samples (no p-values)
    sub_se <- processed_se[, SummarizedExperiment::colData(processed_se)$sample_name %in% c("control_01", "hfref_patient_01")]
    deSp_se <- suppressWarnings(
        deSp_twoGroup(
            sub_se, ref_group='ctrl', test='t-test',
            significant='padj', p_cutoff=0.05, FC_cutoff=1, transform='log10')
    )
    res <- plot_deSp_twoGroup(deSp_se)
    expect_deSp_plots(deSp_se, res)
    expect_true(is.character(res$interactive_volcanoPlot))
    expect_true(is.character(res$static_volcanoPlot))
})

test_that("plot_deSp_twoGroup handles errors correctly", {
    # Test with invalid input
    expect_error(plot_deSp_twoGroup("not a SummarizedExperiment object"),
                 "The input must be a SummarizedExperiment object construct by upstream analysis function, e.g. deSp_twoGroup, deChar_twoGroup...")
})
