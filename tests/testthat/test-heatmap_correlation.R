library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("profiling_data")
processed_se <- data_process(
    profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage')

# get char list
char_list <- list_lipid_char(processed_se)$common_list
names(char_list) <- NULL

expect_corr_heatmap <- function(result) {
    if (!is.character(result)) {
        # Check if `res` is a list
        expect_true(is.list(result), info = "Output should be a list")
        expect_true("interactive_heatmap" %in% names(result), info = "'interactive_heatmap' should be in the list")
        expect_true("static_heatmap" %in% names(result), info = "'static_heatma' should be in the list")
        expect_true("corr_coef_matrix" %in% names(result), info = "'corr_coef_matrix' should be in the list")
        expect_s4_class(result$interactive_heatmap, "IheatmapHorizontal")
        expect_s3_class(result$static_heatmap, "recordedplot")
        expect_true(is.matrix(result$corr_coef_matrix))
    }
}

## corr_method: pearson, spearman
## distfun: pearson, spearman, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski
## hclustfun: complete, single, median, average, ward.D, ward.D2, mcquitty, median, centroid

# Test for basic functionality and output structure
test_that("heatmap_correlation returns correct heatmap result", {
    result_sample <- heatmap_correlation(
        processed_se, char=NULL, transform='log10', correlation='pearson',
        distfun='maximum', hclustfun='average', type='sample')
    expect_corr_heatmap(result_sample)
    result_class <- heatmap_correlation(
        processed_se, char='class', transform='log10', correlation='pearson',
        distfun='maximum', hclustfun='average', type='class')
    expect_corr_heatmap(result_class)
})

# Test for heatmap_correlation handles all char
test_that("heatmap_correlation handles all char.", {
    char_list <- char_list[!(char_list %in% c("FA.OH", "Bond.type"))]
    for (i in seq(char_list)) {
        result <- heatmap_correlation(
            processed_se, char=char_list[i], transform='log10', correlation='pearson',
            distfun='maximum', hclustfun='average', type='class')
        expect_corr_heatmap(result)
    }
    char_list <- c("FA.OH", "Bond.type")
    for (i in seq(char_list)) {
        expect_error(heatmap_correlation(
            processed_se, char=char_list[i], transform='log10', correlation='pearson',
            distfun='maximum', hclustfun='average', type='class'),
            "Not enough data for plotting the heatmap.")
    }
})

# Test for heatmap_correlation handles different correlation methods
test_that("heatmap_correlation handles different correlation methods.", {
    result_pearson <- heatmap_correlation(
        processed_se, char=NULL, transform='log10', correlation='pearson',
        distfun='maximum', hclustfun='average', type='sample')
    result_spearman <- heatmap_correlation(
        processed_se, char=NULL, transform='log10', correlation='spearman',
        distfun='maximum', hclustfun='average', type='sample')
    expect_false(identical(result_pearson$corr_coef_matrix, result_spearman$corr_coef_matrix))
})

# Test all combinations of distfun and hclustfun
test_that("heatmap_correlation handles all the distfun & hclustfun  combinations", {
    distfun_list <-  c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
    hclustfun_list <- c("complete", "single", "median", "average", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
    for (distfun in distfun_list) {
        for (hclustfun in hclustfun_list) {
            result <- heatmap_correlation(
                processed_se, char=NULL, transform='log10', correlation='pearson',
                distfun=distfun, hclustfun=hclustfun, type='sample')
            expect_corr_heatmap(result)
        }
    }
})

# Test error handling for invalid parameters
test_that("heatmap_correlation handles invalid parameters correctly", {
    expect_error(heatmap_correlation(
        processed_se, char=NULL, transform='log10', correlation='pearson',
        distfun='maximum', hclustfun='average', type='class'),
        "For 'class' type, char cannot be NULL.")
    expect_error(
        heatmap_correlation(
            processed_se, char=12345, transform='log10', correlation='pearson',
            distfun='maximum', hclustfun='average', type='class'),
        "Wrong char input, you can view the available char list by list_lipid_char function.")
    expect_error(
        heatmap_correlation(
            processed_se, char=NULL, transform=NULL, correlation='pearson',
            distfun='maximum', hclustfun='average', type='sample'),
        "transform must be one of 'none', 'log10', 'square', or 'cube'.")
    expect_error(
        heatmap_correlation(
            processed_se, char=NULL, transform='log10', correlation="ABC",
            distfun='maximum', hclustfun='average', type='sample'),
        "correlation must be one of 'pearson' or 'spearman'.")
    expect_error(
        heatmap_correlation(
            processed_se, char=NULL, transform='log10', correlation="pearson",
            distfun='ward.D2', hclustfun='average', type='sample'),
        "distfun must be one of 'pearson', 'spearman', 'kendall', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', or 'minkowski'.")
    expect_error(
        heatmap_correlation(
            processed_se, char=NULL, transform='log10', correlation="pearson",
            distfun='maximum', hclustfun='class', type='sample'),
        "hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    expect_error(
        heatmap_correlation(
            processed_se, char=NULL, transform='log10', correlation="pearson",
            distfun='maximum', hclustfun='average', type='lipid'),
        "type must be one of 'sample' or 'class'.")
})

# Test handling of unprocessed data
test_that("heatmap_correlation handles unprocessed data correctly", {
    expect_error(
        heatmap_correlation(
            profiling_data, char=NULL, transform='log10',
            correlation="pearson", distfun='maximum', hclustfun='average',
            type='sample'),
        "Detect 0 values or NAs in abundance data. Please perform data imputation by the data_process function first."
    )
})

# Test different transform methods
test_that("heatmap_correlation handles different transform methods", {
    transform_methods <- c('none', 'log10', 'square', 'cube')
    for (method in transform_methods) {
        result <- heatmap_correlation(
            processed_se, char=NULL, transform=method,
            correlation='pearson', distfun='maximum', hclustfun='average',
            type='sample')
        expect_corr_heatmap(result)
    }
})

