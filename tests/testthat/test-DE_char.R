## Load data and create example data sets
data(DE_data)

## Start tests
test_that("DE_char function works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    ## get lipid characteristics
    lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(DE_data))
    char_var <- colnames(lipid_char_table)[-1]
    ## aggregated(sum) expression data by selected characteristics
    Spe2Char_result <- Species2Char(DE_data, char_var = char_var[4])
    ## data processing of exp_data (without log10 transformation)
    exp_transform_non_log <- data_process(
        Spe2Char_result, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
        NA2what='min', ymin=0.5, pct_transform=TRUE, data_transform=FALSE,
        trans_type='log', centering=FALSE, scaling=FALSE)
    ## conduct deferentially expressed of lipid characters
    DE_char_result <- DE_char(
        exp_transform_non_log, data_transform=TRUE, paired=FALSE,
        sig_pvalue=0.05, sig_FC=2, insert_ref_group=NULL, ref_group=NULL)
    ## test output
    expect_s4_class(DE_char_result$char_exp_data, "SummarizedExperiment")
    expect_s4_class(DE_char_result$char_table_all, "SummarizedExperiment")
    expect_s4_class(DE_char_result$combined_table, "SummarizedExperiment")
    expect_s4_class(
        DE_char_result$combine_result_table, "SummarizedExperiment")
    expect_s4_class(DE_char_result$bar_table, "SummarizedExperiment")
    expect_s4_class(DE_char_result$char_box, "SummarizedExperiment")
    expect_s4_class(DE_char_result$char_table_all_sig, "SummarizedExperiment")
    ## test paired=TRUE
    DE_char_result <- DE_char(
        exp_transform_non_log, data_transform=TRUE, paired=TRUE,
        sig_pvalue=0.05, sig_FC=2, insert_ref_group=NULL, ref_group=NULL)
})
