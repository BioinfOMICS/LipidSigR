## Load data and create example data sets
data(profiling_data)

## Start tests
test_that("corr_heatmap works", {
    expect_s4_class(profiling_data, "SummarizedExperiment")
    ## data processing of exp_data
    exp_transform_SE <- data_process(
        profiling_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
        NA2what='min', ymin=0.5,  pct_transform=TRUE, data_transform=FALSE,
        trans_type='log', centering=FALSE, scaling=FALSE)
    ## sample correlation calculation
    corr_sample_sample <- corr_heatmap(
        exp_transform_SE, corr_method="pearson", distfun="maximum",
        hclustfun="average", type='sample')
    corr_sample_lipid <- corr_heatmap(
        exp_transform_SE, corr_method="pearson", distfun="maximum",
        hclustfun="average", type='lipid')
    ## test output
    expect_s4_class(corr_sample_sample, "SummarizedExperiment")
    expect_s4_class(corr_sample_lipid, "SummarizedExperiment")
    })

