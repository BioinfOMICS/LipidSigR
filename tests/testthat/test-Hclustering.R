# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)
test_groupInfo <- DE_group_info[,
                                c("pair", "sample_name", "group", "label_name")]

# Start tests
test_that("Hclustering function works", {
  exp_transform_table <- data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
  exp_transform_non_log <- data_process(DE_exp_data, exclude_var_missing=TRUE,
     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
     data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
  lipid_char_filter <-
    DE_lipid_char_table[DE_lipid_char_table$feature %in%
                          exp_transform_table$feature, ]
  test_lipid_char <- data.frame(feature = lipid_char_filter$feature,
                                class = lipid_char_filter$totallength,
                                totallength = lipid_char_filter$class,
                                totaldb = lipid_char_filter$class,
                                totaloh = lipid_char_filter$class)
  char_var <- colnames(lipid_char_filter)[-1]
  DE_species_table_sig <- DE_species_2(exp_transform_non_log,
      data_transform=TRUE, group_info = DE_group_info, paired=FALSE,
      test='t.test', adjust_p_method='BH', sig_stat='p.adj', sig_pvalue=0.05,
      sig_FC=2, plotly=FALSE)$DE_species_table_sig
  # test exp_data at least 2 samples.
  expect_error(Hclustering(exp_transform_table[, 1], DE_species_table_sig,
      DE_group_info, lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  # test exp_data number of lipids names (features) must be more than 2.
  expect_error(Hclustering(exp_transform_table[1, ], DE_species_table_sig,
      DE_group_info, lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  # test group_info
  expect_error(Hclustering(exp_transform_table[1, ], DE_species_table_sig,
      test_groupInfo, lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  expect_error(Hclustering(exp_transform_table[1, ], DE_species_table_sig,
      test_groupInfo[, 1:3], lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  # test lipid_char_table content of column 'class' must be characters
  expect_error(Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, test_lipid_char, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  test_lipid_char$class <- lipid_char_filter$class
  # test lipid_char_table content of column 'totallength' must be numeric
  expect_error(Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, test_lipid_char, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  test_lipid_char$totallength <- lipid_char_filter$totallength
  # test lipid_char_table content of column 'totaldb' must be numeric
  expect_error(Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, test_lipid_char, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  test_lipid_char$totaldb <- lipid_char_filter$totaldb
  # test lipid_char_table content of column 'totaloh' must be numeric
  expect_error(Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, test_lipid_char, char_var = char_var[1],
      istfun='pearson', hclustfun='complete', plotly=TRUE, type='all'))
  # all
  result_all <- Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='all')
  expect_type(result_all$data, "double")
  expect_s4_class(result_all$heatmap, "Iheatmap")
  result_all_ggplot <- Hclustering(exp_transform_table,
      DE_species_table_sig, DE_group_info, lipid_char_filter,
      char_var = char_var[1], distfun='pearson', hclustfun='complete',
      plotly=FALSE, type='all')
  expect_s3_class(result_all_ggplot$heatmap, "recordedplot")

  # sig
  result_sig <- Hclustering(exp_transform_table, DE_species_table_sig,
      DE_group_info, lipid_char_filter, char_var = char_var[1],
      distfun='pearson', hclustfun='complete', plotly=TRUE, type='sig')
  expect_type(result_sig$data, "double")
  expect_s4_class(result_sig$heatmap, "Iheatmap")
  result_sig_ggplot <- Hclustering(exp_transform_table,
      DE_species_table_sig, DE_group_info, lipid_char_filter,
      char_var = char_var[1], distfun='pearson', hclustfun='complete',
      plotly=FALSE, type='sig')
  expect_s3_class(result_sig_ggplot$heatmap, "recordedplot")
})
