## Load data and create example data sets
data(DE_data)

## Start tests
test_that("DE_sub_char functions works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    ## get lipid characteristics
    lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(DE_data))
    char_var <- colnames(lipid_char_table)[-1]
    ## subgroup deferentially expressed of lipid characters
    DE.sub.char <- DE_sub_char(
        DE_data, data_transform=TRUE, split_var=char_var[2],
        char_var=char_var[4], paired = FALSE, sig_pvalue=0.05, sig_FC=2,
        exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
        zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
        pct_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
    ## test output
    expect_s4_class(DE.sub.char$result_table1, "SummarizedExperiment")
    expect_s4_class(DE.sub.char$result_table2, "SummarizedExperiment")
    expect_s4_class(DE.sub.char$result_table3, "SummarizedExperiment")
    expect_s4_class(DE.sub.char$result_table4, "SummarizedExperiment")
})
