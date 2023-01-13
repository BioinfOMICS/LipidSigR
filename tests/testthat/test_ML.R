library(SHAPforxgboost)
# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)

# Start tests
test_that("input data are data frames.", {
  expect_s3_class(ML_exp_data, "data.frame")
  expect_s3_class(ML_lipid_char_table, "data.frame")
  expect_s3_class(ML_condition_table, "data.frame")
  expect_equal(nrow(ML_exp_data), nrow(ML_lipid_char_table))
  expect_equal(ncol(ML_exp_data)-1, nrow(ML_condition_table))
})

test_that("data_process() function computed successfullly.", {
  expect_s3_class(ML_exp_data, "data.frame")
  exp_transform_table <- data_process(ML_exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=TRUE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
  expect_s3_class(exp_transform_table, "data.frame")
  expect_equal(ncol(exp_transform_table), ncol(ML_exp_data))
})

test_that("ML plot function computed successfullly.", {
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
  ## ROC plot
  ROC_result <- ROC_plot_all(ML_output[[3]], ML_output[[5]],
                             feature_n=10, plotly=TRUE)
  expect_s3_class(ROC_result[[3]], "data.frame")
  expect_s3_class(ROC_result[[1]], "data.frame")
  expect_s3_class(ROC_result[[2]], "plotly")
  expect_s3_class(ROC_result[[4]], "plotly")
  ## PR plot
  PR_result <- PR_plot_all(ML_output[[4]], ML_output[[5]],
                           feature_n=10, plotly=TRUE)
  expect_s3_class(PR_result[[3]], "data.frame")
  expect_s3_class(PR_result[[1]], "data.frame")
  expect_s3_class(PR_result[[2]], "plotly")
  expect_s3_class(PR_result[[4]], "plotly")
  ## model evaluation plot
  evaluate_result <- evalution_plot(ML_output[[2]], method='Accuracy',
                                    plotly=TRUE)
  expect_s3_class(evaluate_result[[1]], "data.frame")
  expect_s3_class(evaluate_result[[2]], "plotly")
  ## probability plot
  prob_result <- probability_plot(ML_output[[1]], feature_n=10,
                                  plotly=TRUE)
  expect_s3_class(prob_result[[1]], "data.frame")
  expect_s3_class(prob_result[[4]], "data.frame")
  expect_s3_class(prob_result[[2]], "plotly")
  expect_s3_class(prob_result[[3]], "plotly")
})

test_that("Feature importance function computed successfullly.", {
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
  ## Algorithm-based
  feature_result <- feature_plot(ML_output[[6]], ML_output[[7]],
                                 feature_n=10, nfold=10,
                                 plotly=TRUE)
  expect_s3_class(feature_result[[1]], "data.frame")
  expect_s3_class(feature_result[[3]], "data.frame")
  expect_s3_class(feature_result[[2]], "plotly")
  expect_s3_class(feature_result[[4]], "plotly")

  ## SHAP
  SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
                      best_model_feature=ML_output[[9]],
                      ML_method='Random_forest', feature_n=10, nsim=5,
                      plotly=TRUE)
  expect_s3_class(SHAP_output[[1]], "data.frame")
  expect_s3_class(SHAP_output[[2]], "data.frame")
  expect_s3_class(SHAP_output[[3]], "plotly")
  expect_s3_class(SHAP_output[[4]], "plotly")

  ## SHAP feature importance of 10 samples
  SHAP_sample_result <- SHAP_sample(SHAP_output[[2]], n_sample=10,
                                    plotly=TRUE)
  expect_s3_class(SHAP_sample_result, "plotly")
  ## SHAP forceplot
  SHAP_force_result <- SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
                                      cluster_method="ward.D", group_num=10,
                                      plotly=TRUE)
  expect_s3_class(SHAP_force_result[[1]], "data.frame")
  expect_s3_class(SHAP_force_result[[2]], "plotly")
  ## SHAP dependence plot
  SHAP_depend_result <- SHAP_dependence_plot(SHAP_output[[2]], x="C38.6.PC",
                                             y="C38.6.PC", color_var="C38.6.PC",
                                             plotly=TRUE)
  expect_s3_class(SHAP_depend_result, "plotly")
})

test_that("ML network function computed successfullly.", {
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
  model_net <- model_for_net(ML_data[[2]], ML_method='Random_forest',
                             varimp_method='Algorithm-based', ML_output[[8]],
                             ML_output[[9]], feature_num=10, nsim=5)
  cor_net_result <- cor_network(ML_data[[1]], ML_lipid_char_table,
                                model_net[[2]], model_net[[3]],
                                cor_method='pearson', edge_cutoff=0)
  expect_s3_class(cor_net_result, "visNetwork")
})
