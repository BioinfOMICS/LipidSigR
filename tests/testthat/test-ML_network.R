# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)

# Start tests
test_that("network functions works", {
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
  expect_error(model_for_net(ML_data[[2]], ML_method='Random_forest',
      varimp_method='Algorithmbased', ML_output[[8]], ML_output[[9]],
      feature_num=10, nsim=5))
  expect_error(model_for_net(ML_data[[2]], ML_method='Randomforest',
      varimp_method='Algorithm-based', ML_output[[8]], ML_output[[9]],
      feature_num=10, nsim=5))
  expect_error(model_for_net(ML_data[[2]], ML_method='Random_forest',
      varimp_method='Algorithm-based', ML_output[[8]], ML_output[[9]],
      feature_num="10", nsim=5))
  expect_error(model_for_net(ML_data[[2]], ML_method='Random_forest',
      varimp_method='Algorithm-based', ML_output[[8]], ML_output[[9]],
      feature_num=10, nsim="5"))
  model_for_net(ML_data[[2]], ML_method='Random_forest',
                varimp_method='SHAP analysis', ML_output[[8]], ML_output[[9]],
                feature_num=10, nsim=5)
  model_net <- model_for_net(ML_data[[2]], ML_method='Random_forest',
      varimp_method='Algorithm-based', ML_output[[8]], ML_output[[9]],
      feature_num=10, nsim=5)

  ## cor_network
  # test exp_transform_table number of lipids names must be more than 10.
  expect_error(cor_network(ML_data[[1]][1:9,], ML_lipid_char_table,
      model_net[[2]], model_net[[3]], cor_method='pearson', edge_cutoff=0,
      plotly=TRUE))
  # test exp_transform_table at least 60 samples.
  expect_error(cor_network(ML_data[[1]][, 1:59], ML_lipid_char_table,
      model_net[[2]], model_net[[3]], cor_method='pearson', edge_cutoff=0,
      plotly=TRUE))
  # test exp_transform_table first column must contain a list of lipids names
  expect_error(cor_network(ML_data[[1]][, -1], ML_lipid_char_table,
      model_net[[2]], model_net[[3]], cor_method='pearson', edge_cutoff=0,
      plotly=TRUE))
  cor_net_result <- cor_network(ML_data[[1]], ML_lipid_char_table,
      model_net[[2]], model_net[[3]], cor_method='pearson', edge_cutoff=0,
      plotly=TRUE)
  expect_s3_class(cor_net_result$visNetwork_plot, "visNetwork")
  cor_net_ggplot <- cor_network(ML_data[[1]], ML_lipid_char_table,
      model_net[[2]], model_net[[3]], cor_method='pearson', edge_cutoff=0,
      plotly=FALSE)
  expect_s3_class(cor_net_ggplot$visNetwork_plot, "ggplot")
})
