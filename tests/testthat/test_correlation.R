# Load data and create example data sets
data(corr_exp_data)
data(corr_lipid_char_table)
data(corr_condition_table)
data(corr_adjusted_table)

# Start tests
test_that("input data are data frames.", {
  expect_s3_class(corr_exp_data, "data.frame")
  expect_s3_class(corr_lipid_char_table, "data.frame")
  expect_s3_class(corr_condition_table, "data.frame")
  expect_s3_class(corr_adjusted_table, "data.frame")
  expect_equal(nrow(corr_exp_data), nrow(corr_lipid_char_table))
  expect_equal(ncol(corr_exp_data)-1, nrow(corr_condition_table))
  expect_equal(ncol(corr_exp_data)-1, nrow(corr_adjusted_table))
})

test_that("data_process() function computed successfullly.", {
  exp_transform_table <- data_process(corr_exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=TRUE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
  expect_s3_class(exp_transform_table, "data.frame")
  expect_equal(ncol(exp_transform_table), ncol(corr_exp_data))
})

test_that("Clin_Cor_heatmap() function computed successfullly.", {
  exp_transform <- data_process(corr_exp_data, exclude_var_missing=TRUE,
                                missing_pct_limit=50,
                                replace_zero=TRUE, zero2what='min',
                                xmin=0.5, replace_NA=TRUE,
                                NA2what='min', ymin=0.5,
                                pct_transform=TRUE, data_transform=TRUE,
                                trans_type='log', centering=FALSE,
                                scaling=FALSE)
  COspec_clinCor <- Clin_Cor_heatmap(exp_transform, corr_condition_table,
                                     test = 'pearson', adjust_p_method = 'BH',
                                     sig_stat = 'p.adj', sig_pvalue=1,
                                     sig_cor_coef=0, heatmap_col='statistic',
                                     distfun='spearman', hclustfun='average',
                                     plotly=TRUE)
  expect_s3_class(COspec_clinCor$Cor_table_all, "data.frame")
  expect_s3_class(COspec_clinCor$Cor_table_sig, "data.frame")
  expect_type(COspec_clinCor$Cor_reorder_mat, "double")
  expect_s3_class(COspec_clinCor$Cor_table_plot, "plotly")
})

test_that("Clin_LR_heatmap() function computed successfullly.", {
  exp_transform <- data_process(corr_exp_data, exclude_var_missing=TRUE,
                                missing_pct_limit=50,
                                replace_zero=TRUE, zero2what='min',
                                xmin=0.5, replace_NA=TRUE,
                                NA2what='min', ymin=0.5,
                                pct_transform=TRUE, data_transform=TRUE,
                                trans_type='log', centering=FALSE,
                                scaling=FALSE)
  COspec_clin_LR <- Clin_LR_heatmap(exp_transform, corr_condition_table,
                                    corr_adjusted_table,
                                    adjust_p_method = 'BH',
                                    sig_stat = 'p.adj', sig_pvalue = 1,
                                    distfun='spearman', hclustfun='centroid',
                                    heatmap_col='beta_coef',
                                    plotly=TRUE)
  expect_s3_class(COspec_clin_LR$LR_table_all, "data.frame")
  expect_s3_class(COspec_clin_LR$LR_table_sig, "data.frame")
  expect_type(COspec_clin_LR$LR_reorder_mat, "double")
  expect_s3_class(COspec_clin_LR$LR_table_plot, "plotly")
})
