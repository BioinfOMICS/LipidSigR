# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)

# Start tests
test_that("Enrichment works", {
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
  expect_error(Enrichment(DE_species_table_sig, test_lipid_char,
      char_var=char_var[1], sig_pvalue=0.05, plotly=TRUE))
  test_lipid_char$class <- lipid_char_filter$class
  # test lipid_char_table content of column 'totallength' must be numeric
  expect_error(Enrichment(DE_species_table_sig, test_lipid_char,
      char_var=char_var[1], sig_pvalue=0.05, plotly=TRUE))
  test_lipid_char$totallength <- lipid_char_filter$totallength
  # test lipid_char_table content of column 'totaldb' must be numeric
  expect_error(Enrichment(DE_species_table_sig, test_lipid_char,
      char_var=char_var[1], sig_pvalue=0.05, plotly=TRUE))
  test_lipid_char$totaldb <- lipid_char_filter$totaldb
  # test lipid_char_table content of column 'totaloh' must be numeric
  expect_error(Enrichment(DE_species_table_sig, test_lipid_char,
      char_var=char_var[1], sig_pvalue=0.05, plotly=TRUE))
  enrich_result <- Enrichment(DE_species_table_sig,
      lipid_char_table = lipid_char_filter, char_var=char_var[1],
      sig_pvalue=0.05, plotly=TRUE)
  expect_s3_class(enrich_result$enrich_char_table, "data.frame")
  expect_s3_class(enrich_result$enrich_char_barplot, "plotly")
  enrich_ggplot <- Enrichment(DE_species_table_sig,
      lipid_char_table = lipid_char_filter, char_var=char_var[1],
      sig_pvalue=0.05, plotly=FALSE)
  expect_s3_class(enrich_ggplot$enrich_char_barplot, "ggplot")
})
