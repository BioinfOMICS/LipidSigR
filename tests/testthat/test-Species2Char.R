## Load data and create example data sets
data(DE_data)
lipid_char_table <- as.data.frame(
    SummarizedExperiment::rowData(DE_data))
char_var <- colnames(lipid_char_table)[-1]

## Start tests
test_that("Species2Char function works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    ## conduct Species2Char
    exp_data_Spe2Char <- Species2Char(DE_data, char_var=char_var[1])
    ## test output
    expect_s4_class(exp_data_Spe2Char, "SummarizedExperiment")
    ## test char of FA
    exp_data_Spe2Char <- Species2Char(DE_data, char_var=char_var[8])
    expect_s4_class(exp_data_Spe2Char, "SummarizedExperiment")
    })
