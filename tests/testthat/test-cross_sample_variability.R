library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

data("profiling_data")

expect_profiling_plots <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("interactive_lipid_number_barPlot" %in% names(res), info = "'interactive_lipid_number_barPlot' should be in the list")
    expect_true("interactive_lipid_amount_barPlot" %in% names(res), info = "'interactive_lipid_amount_barPlot' should be in the list")
    expect_true("interactive_lipid_distribution" %in% names(res), info = "'interactive_lipid_distribution' should be in the list")
    expect_true("static_lipid_number_barPlot" %in% names(res), info = "'static_lipid_number_barPlot' should be in the list")
    expect_true("static_lipid_amount_barPlot" %in% names(res), info = "'static_lipid_amount_barPlot' should be in the list")
    expect_true("static_lipid_distribution" %in% names(res), info = "'static_lipid_distribution' should be in the list")
    expect_true("table_total_lipid" %in% names(res), info = "'table_total_lipid' should be in the list")

    # Check char_var & print_static parameter
    expect_s3_class(res$interactive_lipid_number_barPlot, "plotly")
    expect_s3_class(res$interactive_lipid_amount_barPlot, "plotly")
    expect_s3_class(res$interactive_lipid_distribution, "plotly")
    expect_s3_class(res$static_lipid_number_barPlot, "ggplot")
    expect_s3_class(res$static_lipid_amount_barPlot, "ggplot")
    expect_s3_class(res$static_lipid_distribution, "ggplot")
    expect_s3_class(res$table_total_lipid, "data.frame")
}

# Test for basic functionality and output structure
test_that("cross_sample_variability returns a SummarizedExperiment object with expected metadata", {
    result <- cross_sample_variability(profiling_data)
    expect_profiling_plots(result)
})


# test_that("cross_sample_variability can handle special condition", {
#     # negative value
#     sub <- profiling_data[1:3, ]
#     # result <- cross_sample_variability(sub)
#     # expect_profiling_plots(result)
#     sub <- sub[, SummarizedExperiment::colData(sub)$sample_name %in% c("control_01", "hfref_patient_01")]
#     result <- cross_sample_variability(sub)
#     expect_profiling_plots(result)
# })
