# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)

# Start tests
test_that("multiplication works", {
  char_var <- colnames(DE_lipid_char_table)[-1]
  exp_data_Spe2Char <- Species2Char(DE_exp_data, DE_lipid_char_table,
                                    char_var = char_var[4])
  exp_transform_non_log <- data_process(exp_data_Spe2Char,
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, data_transform=FALSE, trans_type='log',
      centering=FALSE, scaling=FALSE)
  # test exp_data at least 2 samples.
  expect_error(DE_char_2(exp_transform_non_log[, 1], data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, sig_pvalue=0.05, sig_FC=2,
      insert_ref_group=NULL, ref_group=NULL, plotly=TRUE))
  # data_transform=FALSE
  DE_char_2(exp_transform_non_log, data_transform=FALSE,
      group_info = DE_group_info, paired=FALSE, sig_pvalue=0.05, sig_FC=2,
      insert_ref_group=NULL, ref_group=NULL, plotly=TRUE)

  DE_char_result <- DE_char_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, sig_pvalue=0.05, sig_FC=1,
      insert_ref_group=NULL, ref_group=NULL, plotly=TRUE)
  expect_s3_class(DE_char_result$DE_char_exp_data, "data.frame")
  expect_s3_class(DE_char_result$DE_char_table_all, "data.frame")
  expect_s3_class(DE_char_result$DE_char_combined_table, "data.frame")
  expect_s3_class(DE_char_result$DE_char_combine_result_table, "data.frame")
  expect_s3_class(DE_char_result$DE_char_barplot, "plotly")
  expect_s3_class(DE_char_result$DE_char_barplot_sqrt, "plotly")
  expect_s3_class(DE_char_result$DE_char_trendplot, "plotly")
  expect_s3_class(DE_char_result$DE_char_trendplot_sqrt, "plotly")
  expect_s3_class(DE_char_result$DE_char_boxplot, "plotly")
  DE_char_ggplot <- DE_char_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, sig_pvalue=0.05, sig_FC=1,
      insert_ref_group=NULL, ref_group=NULL, plotly=FALSE)
  expect_s3_class(DE_char_ggplot$DE_char_exp_data, "data.frame")
  expect_s3_class(DE_char_ggplot$DE_char_table_all, "data.frame")
  expect_s3_class(DE_char_ggplot$DE_char_combined_table, "data.frame")
  expect_s3_class(DE_char_ggplot$DE_char_combine_result_table, "data.frame")
  expect_s3_class(DE_char_ggplot$DE_char_barplot, "ggplot")
  expect_s3_class(DE_char_ggplot$DE_char_barplot_sqrt, "ggplot")
  expect_s3_class(DE_char_ggplot$DE_char_trendplot, "ggplot")
  expect_s3_class(DE_char_ggplot$DE_char_trendplot_sqrt, "ggplot")
  expect_s3_class(DE_char_ggplot$DE_char_boxplot, "ggplot")

})
