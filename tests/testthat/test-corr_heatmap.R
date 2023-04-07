# Load data and create example data sets
data(profiling_exp_data)
data(profiling_lipid_char_table)

# Start tests
test_that("corr_heatmap works", {
  exp_transform_table <- data_process(profiling_exp_data,
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, data_transform=TRUE, trans_type='log',
      centering=FALSE, scaling=FALSE)

  # sample
  sample_result <- corr_heatmap(exp_transform_table, corr_method="pearson",
      distfun="maximum", hclustfun="average", plotly=TRUE, type='sample')
  expect_type(sample_result$corr_coef, "double")
  expect_type(sample_result$corr_p, "double")
  expect_type(sample_result$reorder_data, "double")
  expect_s4_class(sample_result$heatmap, "Iheatmap")
  sample_ggplot <- corr_heatmap(exp_transform_table, corr_method="pearson",
      distfun="maximum", hclustfun="average", plotly=FALSE, type='sample')
  expect_s3_class(sample_ggplot$heatmap, "recordedplot")
  # lipid
  lipid_result <- corr_heatmap(exp_transform_table, corr_method="pearson",
      distfun="maximum", hclustfun="average", plotly=TRUE, type='lipid')
  expect_type(lipid_result$corr_coef, "double")
  expect_type(lipid_result$corr_p, "double")
  expect_type(lipid_result$reorder_data, "double")
  expect_s4_class(lipid_result$heatmap, "Iheatmap")
  lipid_ggplot <- corr_heatmap(exp_transform_table, corr_method="pearson",
      distfun="maximum", hclustfun="average", plotly=FALSE, type='lipid')
  expect_s3_class(lipid_ggplot$heatmap, "recordedplot")


})
