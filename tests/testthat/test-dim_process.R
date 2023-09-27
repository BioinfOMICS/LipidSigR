## Load data and create example data sets
data(DE_data)

## Start tests
test_that("dim_process works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    ## data processing of exp_data
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
        NA2what='min', ymin=0.5,  pct_transform=TRUE, data_transform=FALSE,
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
    sig_feature <- as.data.frame(
        SummarizedExperiment::assay(DE_species_result$exp_data_stat_sig))[[1]]
    ## conduct dimension reduction data processing
    dim_process_PCA <- dim_process(
        exp_transform_SE, sig_feature=sig_feature, type='PCA',
        insert_ref_group=NULL, ref_group=NULL)
    dim_process_PLSDA <- dim_process(
        exp_transform_SE, sig_feature=sig_feature, type='PLSDA',
        insert_ref_group=NULL, ref_group=NULL)
    dim_process_result <- dim_process(
        exp_transform_SE, sig_feature=NULL, type='PLSDA',
        insert_ref_group=NULL, ref_group=NULL)

    ## test output
    expect_s4_class(dim_process_PCA, "SummarizedExperiment")
    expect_s4_class(dim_process_PLSDA, "SummarizedExperiment")
    expect_s4_class(dim_process_result, "SummarizedExperiment")
    })
