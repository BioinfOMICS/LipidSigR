# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
test_lipid_char <- data.frame(feature = DE_lipid_char_table$feature,
                              class = DE_lipid_char_table$totallength,
                              totallength = DE_lipid_char_table$class,
                              totaldb = DE_lipid_char_table$class,
                              totaloh = DE_lipid_char_table$class)

# Start tests
test_that("data_process works", {
  expect_s3_class(DE_exp_data, "data.frame")
  # test At least 2 samples.
  expect_error(data_process(DE_exp_data[, 1], exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test is.numeric(zero2what)
  data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what=0.5, xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  # test zero2what=='NA'
  data_process(DE_exp_data, exclude_var_missing=TRUE,
    missing_pct_limit=50, replace_zero=TRUE, zero2what='NA', xmin=0.5,
    replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
    data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  # test NA2what == 'mean'
  data_process(DE_exp_data, exclude_var_missing=TRUE,
    missing_pct_limit=50, replace_zero=TRUE, zero2what='NA', xmin=0.5,
    replace_NA=TRUE, NA2what='mean', ymin=0.5, pct_transform=TRUE,
    data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  # test NA2what == 'median'
  data_process(DE_exp_data, exclude_var_missing=TRUE,
    missing_pct_limit=50, replace_zero=TRUE, zero2what='NA', xmin=0.5,
    replace_NA=TRUE, NA2what='median', ymin=0.5, pct_transform=TRUE,
    data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  # test is.numeric(NA2what)
  data_process(DE_exp_data, exclude_var_missing=TRUE,
    missing_pct_limit=50, replace_zero=TRUE, zero2what='NA', xmin=0.5,
    replace_NA=TRUE, NA2what=0.5, ymin=0.5, pct_transform=TRUE,
    data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)

  exp_transform_table <- data_process(DE_exp_data,
      exclude_var_missing=TRUE, missing_pct_limit=50,
      replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
      NA2what='min', ymin=0.5, pct_transform=TRUE, data_transform=TRUE,
      trans_type='log', centering=TRUE, scaling=TRUE)
  expect_s3_class(exp_transform_table, "data.frame")
  expect_equal(ncol(exp_transform_table), ncol(DE_exp_data))
})
