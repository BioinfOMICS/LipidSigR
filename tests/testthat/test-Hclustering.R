## Load data and create example data sets
data(DE_data)

## Start tests
test_that("Hclustering function works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    ## data processing of exp_data
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
        NA2what='min', ymin=0.5, pct_transform=TRUE, data_transform=TRUE,
        trans_type='log', centering=FALSE, scaling=FALSE)
    ## data processing of exp_data (without log10 transformation)
    exp_transform_non_log <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
        NA2what='min', ymin=0.5, pct_transform=TRUE, data_transform=FALSE,
        trans_type='log', centering=FALSE, scaling=FALSE)
    ## filter significant lipid
    DE_species_result <- DE_species(
        exp_transform_non_log, data_transform=TRUE, paired=FALSE,
        test='t.test', adjust_p_method='BH', sig_stat='p.adj',
        sig_pvalue=0.05, sig_FC=2)
    ## get lipid characteristics
    lipid_char_table <- as.data.frame(
        SummarizedExperiment::rowData(exp_transform_SE))
    char_var <- colnames(lipid_char_table)[-1]
    ## conduct hierarchical clustering of all lipid species
    DE_result_sig <- DE_species_result$exp_data_stat_sig
    Hclustering_all <- Hclustering(
        exp_data_SE=exp_transform_SE, DE_result_sig=DE_result_sig, type='all')
    Hclustering_sig <- Hclustering(
        exp_data_SE=exp_transform_SE, DE_result_sig=DE_result_sig, type='sig')
    ## test output
    expect_s4_class(Hclustering_all, "SummarizedExperiment")
    expect_s4_class(Hclustering_sig, "SummarizedExperiment")
})
