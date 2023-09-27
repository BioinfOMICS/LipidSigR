## Load data and create example data sets
data(profiling_data)
lipid_char_table <- as.data.frame(
    SummarizedExperiment::rowData(profiling_data))
char_var <- colnames(lipid_char_table)[-1]

## Start tests
test_that("lipid_char_table_gather function works.", {
    expect_s4_class(profiling_data, "SummarizedExperiment")
    result_table <- lipid_char_table_gather(
        profiling_data,char_var=char_var[8])
    expect_s4_class(result_table, "SummarizedExperiment")
})
