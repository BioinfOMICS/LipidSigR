# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)

# Start tests
test_that("input data are data frames.", {
  expect_s3_class(DE_exp_data, "data.frame")
  expect_s3_class(DE_lipid_char_table, "data.frame")
  expect_s3_class(DE_group_info, "data.frame")
  expect_equal(nrow(DE_exp_data), nrow(DE_lipid_char_table))
})

test_that("DE species function computed successfullly.", {
  expect_s3_class(DE_exp_data, "data.frame")
  exp_transform_table <- data_process(DE_exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=TRUE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
  expect_s3_class(exp_transform_table, "data.frame")
  expect_equal(ncol(exp_transform_table), ncol(DE_exp_data))
  ## DE_species_2
  exp_transform_non_log <- data_process(DE_exp_data, exclude_var_missing=TRUE,
                                        missing_pct_limit=50, replace_zero=TRUE,
                                        zero2what='min', xmin=0.5,
                                        replace_NA=TRUE, NA2what='min',
                                        ymin=0.5, pct_transform=TRUE,
                                        data_transform=FALSE, trans_type='log',
                                        centering=FALSE, scaling=FALSE)
  DE_species_result <- DE_species_2(exp_transform_non_log,
                                    data_transform=TRUE,
                                    group_info = DE_group_info,
                                    paired=FALSE, test='t.test',
                                    adjust_p_method='BH', sig_stat='p.adj',
                                    sig_pvalue=0.05, sig_FC=2, plotly=TRUE)
  expect_s3_class(DE_species_result$DE_species_table_all, "data.frame")
  expect_s3_class(DE_species_result$DE_species_table_sig, "data.frame")
  expect_s3_class(DE_species_result$DE_species_dotchart_sig, "plotly")
  expect_s3_class(DE_species_result$DE_species_maplot, "plotly")
  expect_s3_class(DE_species_result$DE_species_volcano, "plotly")
  ## PLSDA
  DE_species_table_sig <- DE_species_result$DE_species_table_sig
  DEspec_PLSDA <- PLSDA(exp_transform_table, scaling=TRUE,
                        group_info=DE_group_info, ncomp=2,
                        sig_feature=DE_species_table_sig$feature,
                        cluster_method='group_info',
                        group_num=NULL, var1=NULL, var2=NULL,
                        insert_ref_group=NULL, ref_group=NULL, plotly=TRUE)
  expect_type(DEspec_PLSDA[[1]], "list")
  expect_type(DEspec_PLSDA[[2]], "list")
  expect_s3_class(DEspec_PLSDA[[3]], "plotly")
  expect_s3_class(DEspec_PLSDA[[4]], "plotly")
  ## Hclustering
  lipid_char_filter <-
    DE_lipid_char_table[DE_lipid_char_table$feature %in%
                          exp_transform_table$feature, ]
  char_var <- colnames(lipid_char_filter)[-1]
  DEspec_Hcluster<- Hclustering(exp_transform_table, DE_species_table_sig,
                                DE_group_info, lipid_char_filter,
                                char_var = char_var[1],
                                distfun='pearson', hclustfun='complete',
                                plotly=TRUE)
  expect_type(DEspec_Hcluster$all.lipid.data, "double")
  expect_type(DEspec_Hcluster$sig.lipid.data, "double")
  expect_s4_class(DEspec_Hcluster$all.lipid, "Iheatmap")
  expect_s4_class(DEspec_Hcluster$sig.lipid, "Iheatmap")
  ## Sig lipid feature
  sig_feature_result <- Sig_lipid_feature(DE_species_table_sig,
                                          lipid_char_filter,
                                          char_var[1], sig_FC=2,
                                          plotly=TRUE)
  expect_s3_class(sig_feature_result$barPlot, "plotly")
  expect_s3_class(sig_feature_result$lolipop, "plotly")
  expect_s3_class(sig_feature_result$word, "hwordcloud")
  ## enrichment
  enrich_result <- Enrichment(DE_species_table_sig,
                              lipid_char_table = lipid_char_filter,
                              char_var=char_var[1], sig_pvalue=0.05,
                              plotly=TRUE)
  expect_s3_class(enrich_result$enrich_char_table, "data.frame")
  expect_s3_class(enrich_result$enrich_char_barplot, "plotly")
})

test_that("DE characteristics function", {
  char_var <- colnames(DE_lipid_char_table)[-1]
  expect_type(char_var, "character")
  exp_data_Spe2Char <- Species2Char(DE_exp_data, DE_lipid_char_table,
                                    char_var = char_var[4])
  expect_s3_class(exp_data_Spe2Char, "data.frame")
  expect_equal(ncol(exp_data_Spe2Char), ncol(DE_exp_data))
  exp_transform_non_log <- data_process(exp_data_Spe2Char,
                                        exclude_var_missing=TRUE,
                                        missing_pct_limit=50,
                                        replace_zero=TRUE,
                                        zero2what='min', xmin=0.5,
                                        replace_NA=TRUE, NA2what='min',
                                        ymin=0.5, pct_transform=TRUE,
                                        data_transform=FALSE,
                                        trans_type='log',
                                        centering=FALSE, scaling=FALSE)
  DE_char_result <- DE_char_2(exp_transform_non_log, data_transform=TRUE,
                              group_info = DE_group_info, paired=FALSE,
                              sig_pvalue=0.05, sig_FC=2,
                              insert_ref_group=NULL, ref_group=NULL,
                              plotly=TRUE)
  expect_s3_class(DE_char_result$DE_char_exp_data, "data.frame")
  expect_s3_class(DE_char_result$DE_char_table_all, "data.frame")
  expect_s3_class(DE_char_result$DE_char_combined_table, "data.frame")
  expect_s3_class(DE_char_result$DE_char_combine_result_table, "data.frame")
  expect_s3_class(DE_char_result$DE_char_barplot, "plotly")
  expect_s3_class(DE_char_result$DE_char_barplot_sqrt, "plotly")
  expect_s3_class(DE_char_result$DE_char_trendplot, "plotly")
  expect_s3_class(DE_char_result$DE_char_trendplot_sqrt, "plotly")
  expect_s3_class(DE_char_result$DE_char_boxplot, "plotly")

  ## Subgroup of characteristics
  DE.sub.char.2 <- DE_sub_char_2(DE_exp_data, data_transform=TRUE,
                                 lipid_char_table=DE_lipid_char_table,
                                 split_var = char_var[1],
                                 char_var = char_var[4],
                                 group_info = DE_group_info,
                                 paired=FALSE, sig_pvalue=0.05,
                                 sig_FC=2, exclude_var_missing=TRUE,
                                 missing_pct_limit=50,
                                 replace_zero=TRUE, zero2what='min',
                                 xmin=0.5, replace_NA=TRUE,
                                 NA2what='min', ymin=0.5,
                                 pct_transform=TRUE, trans_type='log',
                                 centering=FALSE, scaling=FALSE)
  expect_s3_class(DE.sub.char.2[[1]], "data.frame")
  expect_s3_class(DE.sub.char.2[[2]], "data.frame")
  expect_s3_class(DE.sub.char.2[[3]], "data.frame")
  expect_s3_class(DE.sub.char.2[[4]], "data.frame")
  char.class <- unique(DE.sub.char.2[[2]][1])
  sub_char_result <- DE_sub_char_plot_2(DE.sub.char.2[[2]],
                                        DE.sub.char.2[[3]],
                                        group_info=DE_group_info,
                                        char_var=char_var[4],
                                        split_var=char_var[1],
                                        split_class=char.class[5,],
                                        insert_ref_group=NULL, ref_group=NULL,
                                        plotly=TRUE)
  expect_s3_class(sub_char_result[[1]], "plotly")
  expect_s3_class(sub_char_result[[2]], "plotly")
  expect_s3_class(sub_char_result[[3]], "plotly")
  expect_s3_class(sub_char_result[[4]], "plotly")
  expect_s3_class(sub_char_result[[5]], "plotly")
})
