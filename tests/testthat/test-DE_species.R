## Load data and create example data sets
data("DE_data")

## Start tests
test_that("DE_species function works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE,missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5,
        replace_NA=TRUE, NA2what='min', ymin=0.5,  pct_transform=TRUE,
        data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
    exp_transform_non_log <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5,
        replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
        data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
    DE_species_result <- DE_species(
        exp_transform_non_log, data_transform=TRUE, paired=FALSE,
        test='t.test', adjust_p_method='BH', sig_stat='p.adj',
        sig_pvalue=0.05, sig_FC=2)
    ## test output
    expect_s4_class(DE_species_result$exp_data_stat, "SummarizedExperiment")
    expect_s4_class(
        DE_species_result$exp_data_stat_sig, "SummarizedExperiment")
    ## sig_stat == "p"
    DE_species_result <- DE_species(
        exp_transform_non_log, data_transform=TRUE, paired=FALSE,
        test='t.test', adjust_p_method='BH', sig_stat='p',
        sig_pvalue=0.05, sig_FC=2)
    expect_s4_class(DE_species_result$exp_data_stat, "SummarizedExperiment")
    expect_s4_class(
        DE_species_result$exp_data_stat_sig, "SummarizedExperiment")
    ## test data_transform=FALSE & test='t.test'
    DE_species_result <- DE_species(
        exp_transform_SE, data_transform=FALSE, paired=FALSE,
        test='t.test', adjust_p_method='BH', sig_stat='p.adj',
        sig_pvalue=0.05, sig_FC=2)
    expect_s4_class(DE_species_result$exp_data_stat, "SummarizedExperiment")
    expect_s4_class(
        DE_species_result$exp_data_stat_sig, "SummarizedExperiment")
    ## test data_transform=FALSE & test='wilcoxon test'
    DE_species_result <- DE_species(
        exp_transform_non_log, data_transform=FALSE, paired=FALSE,
        test='wilcox.test', adjust_p_method='BH', sig_stat='p.adj',
        sig_pvalue=0.05, sig_FC=2)
    expect_s4_class(DE_species_result$exp_data_stat, "SummarizedExperiment")
    expect_s4_class(
        DE_species_result$exp_data_stat_sig, "SummarizedExperiment")
    ## test data_transform=TRUE & test='wilcoxon test'
    DE_species_result <- DE_species(
        exp_transform_non_log, data_transform=TRUE, paired=FALSE,
        test='wilcox.test', adjust_p_method='BH', sig_stat='p.adj',
        sig_pvalue=0.05, sig_FC=2)
    expect_s4_class(DE_species_result$exp_data_stat, "SummarizedExperiment")
    expect_s4_class(
        DE_species_result$exp_data_stat_sig, "SummarizedExperiment")
    ## test paired=TRUE
    testthat::expect_warning(
        DE_species(
        exp_transform_non_log, data_transform=TRUE, paired=TRUE, test='t.test',
        adjust_p_method='BH', sig_stat='p', sig_pvalue=0.05, sig_FC=2)
        )
})
