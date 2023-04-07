library(SHAPforxgboost)
# Load data and create example data sets
data(ML_exp_data)
data(ML_lipid_char_table)
data(ML_condition_table)

# Start tests
test_that("SHAP function works", {
  char_var <- colnames(ML_lipid_char_table)[-1]
  ML_data <- ML_data_process(ML_exp_data, group_info = ML_condition_table,
      ML_lipid_char_table, char_var[1], exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
                        ML_method='Random_forest', split_prop=0.3, nfold=10)
  ## SHAP
  expect_error(SHAP(ML_data[[2]], best_model=ML_output[[8]],
      best_model_feature=ML_output[[9]], ML_method='Randomforest',
      feature_n=10, nsim=5, plotly=TRUE))
  expect_error(SHAP(ML_data[[2]], best_model=ML_output[[8]],
      best_model_feature=ML_output[[9]], ML_method='Random_forest',
      feature_n="10", nsim=5, plotly=TRUE))
  expect_error(SHAP(ML_data[[2]], best_model=ML_output[[8]],
      best_model_feature=ML_output[[9]], ML_method='Random_forest',
      feature_n=10, nsim="5", plotly=TRUE))
  SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
      best_model_feature=ML_output[[9]], ML_method='Random_forest',
      feature_n=10, nsim=5, plotly=TRUE)
  expect_s3_class(SHAP_output[[1]], "data.frame")
  expect_s3_class(SHAP_output[[2]], "data.frame")
  expect_s3_class(SHAP_output[[3]], "plotly")
  expect_s3_class(SHAP_output[[4]], "plotly")
  output_ggplot <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
      best_model_feature=ML_output[[9]], ML_method='Random_forest',
      feature_n=10, nsim=5, plotly=FALSE)
  expect_s3_class(output_ggplot[[3]], "ggplot")
  expect_s3_class(output_ggplot[[4]], "ggplot")

  ## SHAP feature importance of 10 samples
  SHAP_sample(SHAP_output[[2]], n_sample=12, plotly=TRUE)
  SHAP_sample_result <- SHAP_sample(SHAP_output[[2]], n_sample=10, plotly=TRUE)
  expect_s3_class(SHAP_sample_result, "plotly")
  SHAP_sample_ggplot <- SHAP_sample(SHAP_output[[2]], n_sample=10, plotly=FALSE)
  expect_s3_class(SHAP_sample_ggplot, "ggplot")

  ## SHAP forceplot
  expect_error(SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
      cluster_method="wardD", group_num=10, plotly=TRUE))
  expect_error(SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
      cluster_method="ward.D", group_num="10", plotly=TRUE))
  expect_warning(SHAP_forceplot(SHAP_output[[1]], topN_feature=12,
      cluster_method="ward.D", group_num=10, plotly=TRUE))
  SHAP_force_result <- SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
      cluster_method="ward.D", group_num=10, plotly=TRUE)
  expect_s3_class(SHAP_force_result[[1]], "data.frame")
  expect_s3_class(SHAP_force_result[[2]], "plotly")
  SHAP_force_ggplot <- SHAP_forceplot(SHAP_output[[1]], topN_feature=10,
      cluster_method="ward.D", group_num=10, plotly=FALSE)
  expect_s3_class(SHAP_force_ggplot[[2]], "ggplot")
  ## SHAP dependence plot
  SHAP_depend_result <- SHAP_dependence_plot(SHAP_output[[2]], x="C38.6.PC",
      y="C38.6.PC", color_var="C38.6.PC", plotly=TRUE)
  expect_s3_class(SHAP_depend_result, "plotly")
  SHAP_depend_ggplot <- SHAP_dependence_plot(SHAP_output[[2]], x="C38.6.PC",
      y="C38.6.PC", color_var="C38.6.PC", plotly=FALSE)
  expect_s3_class(SHAP_depend_ggplot, "ggplot")
})

