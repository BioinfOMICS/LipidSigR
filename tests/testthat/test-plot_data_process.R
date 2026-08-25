library(testthat)

data("profiling_data")
processed_se <- data_process(
    profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

#### Function ####
process_plot_test <- function(res){
    # Verify if it is a list.
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$interactive_boxPlot_before, "plotly")
    expect_s3_class(res$interactive_densityPlot_before, "plotly")
    expect_s3_class(res$interactive_boxPlot_after, "plotly")
    expect_s3_class(res$interactive_densityPlot_after, "plotly")
    expect_s3_class(res$static_boxPlot_before, "ggplot")
    expect_s3_class(res$static_densityPlot_before, "ggplot")
    expect_s3_class(res$static_boxPlot_after, "ggplot")
    expect_s3_class(res$static_densityPlot_after, "ggplot")
}

test_that("Test results from plot_data_process.R", {
    plot <- plot_data_process(profiling_data, processed_se)
    process_plot_test(plot)
})
