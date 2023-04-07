# Load data and create example data sets
data(corr_exp_data)
data(corr_lipid_char_table)
data(corr_condition_table)
data(corr_adjusted_table)

# Start tests
test_that("Clin_LR_heatmap works", {
  exp_transform_table <- data_process(corr_exp_data, exclude_var_missing=TRUE,
      missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
      replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
      data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  # test At least 10 samples.
  expect_error(Clin_LR_heatmap(exp_transform_table[, 1], corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=TRUE))
  # test At least 2 conditions.
  expect_error(Clin_LR_heatmap(exp_transform_table, corr_condition_table[, 1],
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=TRUE))
  # test At least 10 samples.
  expect_error(Clin_LR_heatmap(exp_transform_table, corr_condition_table[1, ],
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=TRUE))
  # test distfun
  expect_error(Clin_LR_heatmap(exp_transform_table, corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='abc', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=TRUE))
  # test hclustfun
  expect_error(Clin_LR_heatmap(exp_transform_table, corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='abc',
      heatmap_col='beta_coef', plotly=TRUE))
  # sig_stat == 'p', heatmap_col == 't_statistic'
  Clin_LR_heatmap(exp_transform_table, corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='t_statistic', plotly=TRUE)

  Clin_LR <- Clin_LR_heatmap(exp_transform_table, corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=TRUE)
  expect_s3_class(Clin_LR$LR_table_all, "data.frame")
  expect_s3_class(Clin_LR$LR_table_sig, "data.frame")
  expect_type(Clin_LR$LR_reorder_mat, "double")
  expect_s3_class(Clin_LR$LR_table_plot, "plotly")
  Clin_LR_ggplot <- Clin_LR_heatmap(exp_transform_table, corr_condition_table,
      corr_adjusted_table, adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue = 1, distfun='spearman', hclustfun='centroid',
      heatmap_col='beta_coef', plotly=FALSE)
  expect_s3_class(Clin_LR_ggplot$LR_table_plot, "recordedplot")
})
