# Load data and create example data sets
data(corr_exp_data)
data(corr_lipid_char_table)
data(corr_condition_table)
data(corr_adjusted_table)

# Start tests
test_that("Clin_Cor_heatmap works", {
  exp_transform_table <- data_process(corr_exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=TRUE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
  # test At least 10 samples.
  expect_error(Clin_Cor_heatmap(exp_transform_table[, 1], corr_condition_table,
      test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue=1, sig_cor_coef=0, heatmap_col='statistic',
      distfun='spearman', hclustfun='average', plotly=TRUE))
  # test At least 2 conditions.
  expect_error(Clin_Cor_heatmap(exp_transform_table, corr_condition_table[, 1],
      test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue=1, sig_cor_coef=0, heatmap_col='statistic',
      distfun='spearman', hclustfun='average', plotly=TRUE))
  # test At least 10 samples.
  expect_error(Clin_Cor_heatmap(exp_transform_table, corr_condition_table[1, ],
      test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue=1, sig_cor_coef=0, heatmap_col='statistic',
      distfun='spearman', hclustfun='average', plotly=TRUE))
  # sig_stat == 'p', heatmap_col == 'cor_coef'
  Clin_Cor_heatmap(exp_transform_table, corr_condition_table,
      test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p',
      sig_pvalue=1, sig_cor_coef=0, heatmap_col='cor_coef',
      distfun='spearman', hclustfun='average', plotly=TRUE)

  ClinCor <- Clin_Cor_heatmap(exp_transform_table, corr_condition_table,
      test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p.adj',
      sig_pvalue=1, sig_cor_coef=0, heatmap_col='statistic',
      distfun='spearman', hclustfun='average', plotly=TRUE)
  expect_s3_class(ClinCor$Cor_table_all, "data.frame")
  expect_s3_class(ClinCor$Cor_table_sig, "data.frame")
  expect_type(ClinCor$Cor_reorder_mat, "double")
  expect_s3_class(ClinCor$Cor_table_plot, "plotly")
  ClinCor_ggplot <- Clin_Cor_heatmap(exp_transform_table, corr_condition_table,
     test = 'pearson', adjust_p_method = 'BH', sig_stat = 'p.adj',
     sig_pvalue=1, sig_cor_coef=0, heatmap_col='statistic',
     distfun='spearman', hclustfun='average', plotly=FALSE)
  expect_s3_class(ClinCor_ggplot$Cor_table_plot, "recordedplot")
})
