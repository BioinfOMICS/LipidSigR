# Load data and create example data sets
data(DE_exp_data)
data(DE_lipid_char_table)
data(DE_group_info)
test_groupInfo <- DE_group_info[,
                                c("pair", "sample_name", "group", "label_name")]

# Start tests
test_that("PLSDA works", {
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
  # test At least 6 samples.
  expect_error(PLSDA(exp_transform_table[, 1], scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='group_info',
        group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE))
  # test At least 6 features.
  expect_error(PLSDA(exp_transform_table[1, ], scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='group_info',
        group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE))
  # test At least 2 significant lipid features.
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature[1],
        cluster_method='group_info',
        group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE))
  # test group info
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=test_groupInfo, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='group_info',
        group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE))
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=test_groupInfo[, 1:3], ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='group_info',
        group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE))
  # test cluster_method=='kmeans'
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='kmeans', group_num=11,
        var1=NULL, var2=NULL, insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='kmeans', group_num=NULL,
        var1=NULL, var2=NULL, insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  PLSDA(exp_transform_table, scaling=TRUE, group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='kmeans',
        group_num=2, var1=NULL, var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE)
  # test cluster_method == 'kmedoids'
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='kmedoids', group_num=NULL,
        var1='abc', var2=NULL, insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  PLSDA(exp_transform_table, scaling=TRUE, group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='kmedoids',
        group_num=2, var1='euclidean', var2=NULL, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE)
  # test cluster_method == 'hclustering'
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='hclustering', group_num=NULL,
        var1='abc', var2=NULL, insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='hclustering', group_num=NULL,
        var1='pearson', var2='abc', insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  PLSDA(exp_transform_table, scaling=TRUE, group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='hclustering',
        group_num=2, var1='pearson', var2='ward.D', insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE)
  # test cluster_method == 'dbscan'
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
        group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature,
        cluster_method='dbscan', group_num=NULL,
        var1='abc', var2=NULL, insert_ref_group=NULL, ref_group=NULL,
        plotly=TRUE))
  expect_error(PLSDA(exp_transform_table, scaling=TRUE,
       group_info=DE_group_info, ncomp=2,
       sig_feature=DE_species_table_sig$feature,
       cluster_method='dbscan', group_num=NULL,
       var1=2, var2='abc', insert_ref_group=NULL, ref_group=NULL,
       plotly=TRUE))
  PLSDA(exp_transform_table, scaling=TRUE, group_info=DE_group_info, ncomp=2,
        sig_feature=DE_species_table_sig$feature, cluster_method='dbscan',
        group_num=2, var1=2, var2=2, insert_ref_group=NULL,
        ref_group=NULL, plotly=TRUE)

  PLSDA_result <- PLSDA(exp_transform_table, scaling=TRUE,
      group_info=DE_group_info, ncomp=2,
      sig_feature=DE_species_table_sig$feature, cluster_method='group_info',
      group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
      ref_group=NULL, plotly=TRUE)
  expect_type(PLSDA_result[[1]], "list")
  expect_type(PLSDA_result[[2]], "list")
  expect_s3_class(PLSDA_result[[3]], "plotly")
  expect_s3_class(PLSDA_result[[4]], "plotly")
  PLSDA_ggplot <- PLSDA(exp_transform_table, scaling=TRUE,
      group_info=DE_group_info, ncomp=2,
      sig_feature=DE_species_table_sig$feature, cluster_method='group_info',
      group_num=NULL, var1=NULL, var2=NULL, insert_ref_group=NULL,
      ref_group=NULL, plotly=FALSE)
  expect_s3_class(PLSDA_ggplot[[3]], "ggplot")
  expect_s3_class(PLSDA_ggplot[[4]], "ggplot")
})
