## Load data and create example data sets
data(profiling_data)
lipid_char_table <- as.data.frame(
    SummarizedExperiment::rowData(profiling_data))
char_var <- colnames(lipid_char_table)[-1]

## Start tests
test_that("exp_profiling function works.", {
    expect_s4_class(profiling_data, "SummarizedExperiment")
    ## conduct profiling function
    exp_profiling_result <- exp_profiling(profiling_data, char_var[1])
    ## test output
    expect_s4_class(exp_profiling_result$tot.num.lip, "SummarizedExperiment")
    expect_s4_class(exp_profiling_result$dens.lip, "SummarizedExperiment")
    expect_s4_class(
        exp_profiling_result$exp.compo.by.lipid, "SummarizedExperiment")
})
