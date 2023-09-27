## Load data and create example data sets
data(DE_data)

## Start tests
test_that("data_process works", {
    expect_s4_class(DE_data, "SummarizedExperiment")
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what="min", xmin=0.5,
        replace_NA=TRUE, NA2what="min", ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=FALSE, scaling=FALSE)
    ## test output
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    ## test zero2what='NA'
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='NA', xmin=0.5,
        replace_NA=TRUE, NA2what="min", ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=FALSE, scaling=FALSE)
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    ## test NA2what='mean'
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='NA', xmin=0.5,
        replace_NA=TRUE, NA2what='mean', ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=FALSE, scaling=FALSE)
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    ## test NA2what='median'
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='NA', xmin=0.5,
        replace_NA=TRUE, NA2what='median', ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=FALSE, scaling=FALSE)
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    ## test centering=TRUE
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='NA', xmin=0.5,
        replace_NA=TRUE, NA2what='median', ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=TRUE, scaling=FALSE)
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    ## test scaling=TRUE
    exp_transform_SE <- data_process(
        DE_data, exclude_var_missing=TRUE, missing_pct_limit=50,
        replace_zero=TRUE, zero2what='NA', xmin=0.5,
        replace_NA=TRUE, NA2what='median', ymin=0.5,
        pct_transform=TRUE, data_transform=FALSE, trans_type="log",
        centering=FALSE, scaling=TRUE)
    expect_s4_class(exp_transform_SE, "SummarizedExperiment")
    })
