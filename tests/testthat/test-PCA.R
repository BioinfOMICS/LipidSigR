# Load data and create example data sets
data(profiling_exp_data)
data(profiling_lipid_char_table)
data(DE_group_info)
test_groupInfo <- DE_group_info[,
                                c("pair", "sample_name", "group", "label_name")]

# Start tests
test_that("PCA function works", {
  exp_transform_table <- data_process(profiling_exp_data,
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, data_transform=TRUE, trans_type='log',
      centering=FALSE, scaling=FALSE)
  # test At least 6 samples.
  expect_error(PCA(exp_transform_table[, 1:5], group_info = NULL,
      sig_feature = NULL, scaling=TRUE, centering=TRUE, cluster_method='kmeans',
      group_num=2, var1 = NULL, var2 = NULL, insert_ref_group=NULL,
      ref_group=NULL, n_PC=c(1,2), top_n_feature=10, plotly=TRUE))
  # test At least 2 features.
  expect_error(PCA(exp_transform_table[1, ], group_info = NULL,
      sig_feature = NULL, scaling=TRUE, centering=TRUE, cluster_method='kmeans',
      group_num=2, var1 = NULL, var2 = NULL, insert_ref_group=NULL,
      ref_group=NULL, n_PC=c(1,2), top_n_feature=10, plotly=TRUE))
  PCA(exp_transform_table, group_info = DE_group_info, sig_feature = NULL,
      scaling=TRUE, centering=TRUE, cluster_method='kmeans', group_num=2,
      var1 = NULL, var2 = NULL, insert_ref_group=NULL, ref_group=NULL,
      n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  # group info
  expect_error(PCA(exp_transform_table[1, ], group_info = test_groupInfo,
      sig_feature = NULL, scaling=TRUE, centering=TRUE, cluster_method='kmeans',
      group_num=2, var1 = NULL, var2 = NULL, insert_ref_group=NULL,
      ref_group=NULL, n_PC=c(1,2), top_n_feature=10, plotly=TRUE))
  expect_error(PCA(exp_transform_table[1, ], group_info = test_groupInfo[, 1:3],
       sig_feature = NULL, scaling=TRUE, centering=TRUE,cluster_method='kmeans',
       group_num=2, var1 = NULL, var2 = NULL, insert_ref_group=NULL,
       ref_group=NULL, n_PC=c(1,2), top_n_feature=10, plotly=TRUE))

  # cluster = kmedoids
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='kmedoids', group_num=2, var1 = "maximum", var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  PCA(exp_transform_table, group_info = DE_group_info, sig_feature = NULL,
      scaling=TRUE, centering=TRUE, cluster_method='kmedoids', group_num=2,
      var1 = 'euclidean', var2 = NULL, insert_ref_group=NULL, ref_group=NULL,
      n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  # cluster = hclustering
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='hclustering', group_num=2, var1 ="abc", var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='hclustering', group_num=2, var1 = 'pearson', var2 = "abc",
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  PCA(exp_transform_table, group_info = DE_group_info, sig_feature = NULL,
      scaling=TRUE, centering=TRUE, cluster_method='hclustering', group_num=2,
      var1 = 'pearson', var2 = "ward.D", insert_ref_group=NULL, ref_group=NULL,
      n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  # cluster = dbscan
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='dbscan', group_num=2, var1 = NULL, var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='dbscan', group_num=2, var1 = 2, var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  PCA(exp_transform_table, group_info = DE_group_info, sig_feature = NULL,
      scaling=TRUE, centering=TRUE, cluster_method='dbscan', group_num=2,
      var1 = 2, var2 = 2, insert_ref_group=NULL, ref_group=NULL,
      n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  # cluster = kmeans
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='kmeans', group_num=0, var1 = NULL, var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))
  expect_error(PCA(exp_transform_table, group_info = DE_group_info,
      sig_feature = NULL, scaling=TRUE, centering=TRUE,
      cluster_method='kmeans', group_num= NULL, var1 = 2, var2 = NULL,
      insert_ref_group=NULL, ref_group=NULL, n_PC=c(1,2), top_n_feature=10,
      plotly=TRUE))

  # test output
  PCA_result <- PCA(exp_transform_table, group_info = NULL, sig_feature = NULL,
      scaling=TRUE, centering=TRUE, cluster_method='kmeans', group_num=2,
      var1 = NULL, var2 = NULL, insert_ref_group=NULL, ref_group=NULL,
      n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  expect_s3_class(PCA_result[[1]], "prcomp")
  expect_s3_class(PCA_result[[2]], "data.frame")
  expect_s3_class(PCA_result[[3]], "data.frame")
  expect_s3_class(PCA_result[[4]], "plotly")
  expect_s3_class(PCA_result[[5]], "plotly")
  expect_s3_class(PCA_result[[6]], "plotly")
  expect_s3_class(PCA_result[[7]], "plotly")
  PCA_ggplot <- PCA(exp_transform_table, group_info = NULL, sig_feature = NULL,
       scaling=TRUE, centering=TRUE, cluster_method='kmeans', group_num=2,
       var1 = NULL, var2 = NULL, insert_ref_group=NULL, ref_group=NULL,
       n_PC=c(1,2), top_n_feature=10, plotly=FALSE)
  expect_s3_class(PCA_ggplot[[4]], "ggplot")
  expect_s3_class(PCA_ggplot[[5]], "ggplot")
  expect_s3_class(PCA_ggplot[[6]], "ggplot")
  expect_s3_class(PCA_ggplot[[7]], "ggplot")
})
