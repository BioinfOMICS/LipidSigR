library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("profiling_data")
processed_se <- data_process(
    profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

# get char_var
char_list <- list_lipid_char(profiling_data)$common_list
names(char_list) <- NULL

expect_profiling_plots <- function(res, char, print_static) {
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("interactive_char_barPlot" %in% names(res), info = "'interactive_char_barPlot' should be in the list")
    expect_true("interactive_lipid_composition" %in% names(res), info = "'interactive_lipid_composition' should be in the list")
    expect_true("static_char_barPlot" %in% names(res), info = "'static_char_barPlot' should be in the list")
    expect_true("static_lipid_composition" %in% names(res), info = "'static_lipid_composition' should be in the list")
    expect_true("table_char_barPlot" %in% names(res), info = "'table_char_barPlot' should be in the list")
    expect_true("table_lipid_composition" %in% names(res), info = "'table_lipid_composition' should be in the list")
    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_char_barPlot, "plotly")
    expect_s3_class(res$interactive_lipid_composition, "plotly")
    expect_s3_class(res$static_char_barPlot, "ggplot")
    expect_s3_class(res$static_lipid_composition, "ggplot")
    expect_s3_class(res$table_char_barPlot, "data.frame")
    expect_s3_class(res$table_lipid_composition, "data.frame")
}

# Test for basic functionality and output structure
test_that("lipid_profiling returns a SummarizedExperiment object with expected metadata", {
    result <- lipid_profiling(processed_se, char_list[1])
    expect_profiling_plots(result)
})

test_that("lipid_profiling function handles all characteristics.", {
    for (i in seq(char_list)) {
        result <- lipid_profiling(processed_se, char=char_list[i])
        expect_profiling_plots(result)
    }
})

test_that("lipid_profiling can handle invalid char", {
    expect_error(
        lipid_profiling(processed_se, char="lipid"),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    expect_error(
        lipid_profiling(processed_se, char=NULL),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
})
