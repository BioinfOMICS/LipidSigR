# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)

# Start tests
test_that("feature_plot funtcion works", {
  char_var <- colnames(ML_lipid_char_table)[-1]
  ML_data <- ML_data_process(ML_exp_data, group_info = ML_condition_table,
                             ML_lipid_char_table, char_var[1],
                             exclude_var_missing=TRUE, missing_pct_limit=50,
                             replace_zero=TRUE, zero2what='min', xmin=0.5,
                             replace_NA=TRUE, NA2what='min', ymin=0.5,
                             pct_transform=TRUE, data_transform=TRUE,
                             trans_type='log', centering=FALSE, scaling=FALSE)
  ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
                        ML_method='Random_forest', split_prop=0.3, nfold=10)
  expect_error(probability_plot(ML_output[[1]], feature_n="10", plotly=TRUE))
  prob_result <- probability_plot(ML_output[[1]], feature_n=10, plotly=TRUE)
  expect_s3_class(prob_result[[1]], "data.frame")
  expect_s3_class(prob_result[[4]], "data.frame")
  expect_s3_class(prob_result[[2]], "plotly")
  expect_s3_class(prob_result[[3]], "plotly")
  result_ggplot <- probability_plot(ML_output[[1]], feature_n=10, plotly=FALSE)
  expect_s3_class(result_ggplot[[2]], "ggplot")
  expect_s3_class(result_ggplot[[3]], "ggplot")
})

