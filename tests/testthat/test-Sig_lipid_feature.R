# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)
test_lipid_char <- data.frame(feature = DE_lipid_char_table$feature,
                              class = DE_lipid_char_table$totallength,
                              totallength = DE_lipid_char_table$class,
                              totaldb = DE_lipid_char_table$class,
                              totaloh = DE_lipid_char_table$class)

# Start tests
test_that("Sig_lipid_feature works", {
  exp_transform_table <- data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  exp_transform_non_log <- data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
  DE_species_result <- DE_species_2(exp_transform_non_log, data_transform=TRUE,
     group_info = DE_group_info, paired=FALSE, test='t.test',
     adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05,
     sig_FC=2, plotly=TRUE)
  DE_species_table_sig <- DE_species_result$DE_species_table_sig
  lipid_char_filter <-
    DE_lipid_char_table[DE_lipid_char_table$feature %in%
                          exp_transform_table$feature, ]
  test_lipid_char <- data.frame(feature = lipid_char_filter$feature,
    class = lipid_char_filter$totallength,
    totallength = lipid_char_filter$class,
    totaldb = lipid_char_filter$class,
    totaloh = lipid_char_filter$class)
  char_var <- colnames(lipid_char_filter)[-1]
  # test lipid_char_table content of column 'class' must be characters
  expect_error(Sig_lipid_feature(DE_species_table_sig, test_lipid_char,
      char_var[1], sig_FC=2, plotly=TRUE))
  test_lipid_char$class <- lipid_char_filter$class
  # test lipid_char_table content of column 'totallength' must be numeric
  expect_error(Sig_lipid_feature(DE_species_table_sig, test_lipid_char,
      char_var[1], sig_FC=2, plotly=TRUE))
  test_lipid_char$totallength <- lipid_char_filter$totallength
  # test lipid_char_table content of column 'totaldb' must be numeric
  expect_error(Sig_lipid_feature(DE_species_table_sig, test_lipid_char,
      char_var[1], sig_FC=2, plotly=TRUE))
  test_lipid_char$totaldb <- lipid_char_filter$totaldb
  # test lipid_char_table content of column 'totaloh' must be numeric
  expect_error(Sig_lipid_feature(DE_species_table_sig, test_lipid_char,
      char_var[1], sig_FC=2, plotly=TRUE))


  sig_feature_result <- Sig_lipid_feature(DE_species_table_sig,
      lipid_char_filter, char_var[1], sig_FC=2, plotly=TRUE)
  expect_s3_class(sig_feature_result$barPlot, "plotly")
  expect_s3_class(sig_feature_result$lolipop, "plotly")
  expect_s3_class(sig_feature_result$word, "hwordcloud")
  sig_feature_ggplot <- Sig_lipid_feature(DE_species_table_sig,
      lipid_char_filter, char_var[1], sig_FC=2, plotly=FALSE)
  expect_s3_class(sig_feature_ggplot$barPlot, "ggplot")
  expect_s3_class(sig_feature_ggplot$lolipop, "ggplot")
  expect_s3_class(sig_feature_ggplot$word, "recordedplot")
})
