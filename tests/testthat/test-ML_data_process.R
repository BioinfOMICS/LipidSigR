# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)
test_lipid_char <- data.frame(feature = ML_lipid_char_table$feature,
                              class = ML_lipid_char_table$totallength,
                              totallength = ML_lipid_char_table$class,
                              totaldb = ML_lipid_char_table$class,
                              totaloh = ML_lipid_char_table$class)

# Start tests
test_that("ML_data_process works", {
  char_var <- colnames(ML_lipid_char_table)[-1]
  # test exp_data first column must contain a list of lipids names
  expect_error(ML_data_process(ML_exp_data[, -1], group_info=ML_condition_table,
     ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test exp_data first column type must be 'character',others must be 'numeric'
  expect_error(ML_data_process(ML_condition_table,group_info=ML_condition_table,
      ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test exp_data at least 60 samples.
  expect_error(ML_data_process(ML_exp_data[, 1:2],group_info=ML_condition_table,
      ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test exp_data number of lipids names (features) must be more than 10.
  expect_error(ML_data_process(ML_exp_data[1:10, ],
      group_info=ML_condition_table, ML_lipid_char_table, char_var[1],
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, data_transform=TRUE, trans_type='log',
      centering=FALSE, scaling=FALSE))
  # test char_var must be included in the lipid_char_table.
  expect_error(ML_data_process(ML_exp_data[, 1:2],group_info=ML_condition_table,
      ML_lipid_char_table, c("a", "b", "c"), exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test lipid_char_table content of column 'class' must be characters
  expect_error(ML_data_process(ML_exp_dat, group_info=ML_condition_table,
      test_lipid_char, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test lipid_char_table content of column 'totallength' must be numeric
  test_lipid_char$class <- ML_lipid_char_table$class
  expect_error(ML_data_process(ML_exp_dat, group_info=ML_condition_table,
      test_lipid_char, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test lipid_char_table content of column 'totaldb' must be numeric
  test_lipid_char$totallength <- ML_lipid_char_table$totallength
  expect_error(ML_data_process(ML_exp_dat, group_info=ML_condition_table,
      test_lipid_char, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))
  # test lipid_char_tablee content of column 'totaloh' must be numeric
  test_lipid_char$totaldb <- ML_lipid_char_table$totaldb
  expect_error(ML_data_process(ML_exp_dat, group_info=ML_condition_table,
      test_lipid_char, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE))

  ML_data <- ML_data_process(ML_exp_data, group_info = ML_condition_table,
      ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)

  expect_s3_class(ML_data[[1]], "data.frame")
  expect_s3_class(ML_data[[2]], "data.frame")


})
