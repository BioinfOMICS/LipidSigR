# Load data and create example data sets
data(DE_exp_data)
data(DE_group_info)
test_exp_data <- data.frame(feature = DE_exp_data$control_01 ,
                            control_01 = DE_exp_data$feature,
                            control_02 = DE_exp_data$control_02,
                            hfref_patient_01 = DE_exp_data$hfref_patient_01,
                            hfref_patient_02 = DE_exp_data$hfref_patient_02)

# Start tests
test_that("DE_species_2 function works", {
  exp_transform_non_log <- data_process(DE_exp_data, exclude_var_missing=TRUE,
                                        missing_pct_limit=50, replace_zero=TRUE,
                                        zero2what='min', xmin=0.5,
                                        replace_NA=TRUE, NA2what='min',
                                        ymin=0.5, pct_transform=TRUE,
                                        data_transform=FALSE, trans_type='log',
                                        centering=FALSE, scaling=FALSE)
  # test exp_data first column must contain a list of lipids names (features)
  expect_error(DE_species_2(test_exp_data, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE))
  # test ncol(exp_data) == 2
  expect_error(DE_species_2(test_exp_data[, 1:2], data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE))
  # test exp_data first column type must be character, others must be numeric
  test_exp_data$feature <- DE_exp_data$feature
  expect_error(DE_species_2(test_exp_data, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE))
  # exp_data at least 2 samples.
  expect_error(DE_species_2(exp_transform_non_log[, 1], data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE))
  # test exp_data number of lipids names (features) must be more than 5
  expect_error(DE_species_2(exp_transform_non_log[1:3, ], data_transform=TRUE,
      group_info = DE_group_info[1:3, ], paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE))
  # test data transform = FALSE & test='t.test'
  DE_species_2(exp_transform_non_log, data_transform=FALSE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE)
  # test data transform = FALSE & test='wilcox.test'
  DE_species_2(exp_transform_non_log, data_transform=FALSE,
      group_info = DE_group_info, paired=FALSE, test='wilcox.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE)
  # test data transform = TRUE & test='wilcox.test'
  DE_species_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='wilcox.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE)
  # sig_stat == "p"
  DE_species_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE)
  DE_species_2(exp_transform_non_log, data_transform=TRUE,
              group_info = DE_group_info, paired=FALSE, test='t.test',
              adjust_p_method='BH', sig_stat='p', sig_pvalue=0.05, sig_FC=2,
              plotly=FALSE)

  DE_species_result <- DE_species_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=TRUE)
  expect_s3_class(DE_species_result$DE_species_table_all, "data.frame")
  expect_s3_class(DE_species_result$DE_species_table_sig, "data.frame")
  expect_s3_class(DE_species_result$DE_species_dotchart_sig, "plotly")
  expect_s3_class(DE_species_result$DE_species_maplot, "plotly")
  expect_s3_class(DE_species_result$DE_species_volcano, "plotly")
  DE_species_ggplot <- DE_species_2(exp_transform_non_log, data_transform=TRUE,
      group_info = DE_group_info, paired=FALSE, test='t.test',
      adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
      plotly=FALSE)
  expect_s3_class(DE_species_ggplot$DE_species_dotchart_sig, "ggplot")
  expect_s3_class(DE_species_ggplot$DE_species_maplot, "ggplot")
  expect_s3_class(DE_species_ggplot$DE_species_volcano, "ggplot")
})
