# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)

# Start tests
test_that("ML model function works", {
  char_var <- colnames(ML_lipid_char_table)[-1]
  ## data process
  ML_data <- ML_data_process(ML_exp_data, group_info = ML_condition_table,
      ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)

  ## ML model
  # test ranking_method
  expect_error(ML_final(ML_data[[2]], ranking_method='ABC',
      ML_method='Random_forest', split_prop=0.3, nfold=3))
  # test ML method
  expect_error(ML_final(ML_data[[2]], ranking_method='Random_forest',
      ML_method='ABC', split_prop=0.3, nfold=3))
  # test nfold must be integer
  expect_error(ML_final(ML_data[[2]], ranking_method='Random_forest',
      ML_method='Random_forest', split_prop=0.3, nfold='3'))
  # test ranking_method == 'p_value'
  ML_final(ML_data[[2]],ranking_method='p_value', ML_method='Random_forest',
           alpha=0.5, split_prop=0.3, nfold=3)
  # test ranking_method == 'pvalue_FC'
  ML_final(ML_data[[2]],ranking_method='pvalue_FC', ML_method='Random_forest',
           split_prop=0.3, nfold=3)
  # test ranking_method == 'Lasso', ML_method == 'Lasso'
  ML_final(ML_data[[2]],ranking_method='Lasso', ML_method='Lasso',
           split_prop=0.3, nfold=3)
  # test ranking_method == 'ElasticNet', ML_method == 'ElasticNet'
  expect_error(ML_final(ML_data[[2]],ranking_method='ElasticNet',
          ML_method='ElasticNet', split_prop=0.3, nfold=3))
  ML_final(ML_data[[2]],ranking_method='ElasticNet', ML_method='ElasticNet',
           alpha=0.5, split_prop=0.3, nfold=3)
  # test ranking_method == 'ROC'
  ML_final(ML_data[[2]],ranking_method='ROC', ML_method='Random_forest',
           split_prop=0.3, nfold=3)
  # test ranking_method == 'SVM', ML_method == 'SVM'
  ML_final(ML_data[[2]],ranking_method='SVM', ML_method='SVM',
           split_prop=0.3, nfold=3)
  # test ranking_method == 'Ridge', ML_method == 'Ridge'
  ML_final(ML_data[[2]],ranking_method='Ridge', ML_method='Ridge',
           split_prop=0.3, nfold=3)
  # test ranking_method == 'Random_forest', ML_method == 'Random_forest'
  ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
      ML_method='Random_forest', split_prop=0.3, nfold=3)
})
