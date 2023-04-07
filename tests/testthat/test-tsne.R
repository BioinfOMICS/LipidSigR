# Load data and create example data sets
data(profiling_exp_data)
data(profiling_lipid_char_table)
data(DE_group_info)
test_groupInfo <- DE_group_info[,
                                c("pair", "sample_name", "group", "label_name")]

# Start tests
test_that("tsne function works", {
  exp_transform_table <- data_process(profiling_exp_data,
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, data_transform=TRUE, trans_type='log',
      centering=FALSE, scaling=FALSE)
  tsne(exp_transform_table, group_info = DE_group_info,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE)
  # test At least 2 samples.
  expect_error(tsne(exp_transform_table[, 1:3], group_info = DE_group_info,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  # test At least 4 features.
  expect_error(tsne(exp_transform_table[1:3, ], group_info = DE_group_info,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  # group info
  expect_error(tsne(exp_transform_table, group_info = test_groupInfo,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  expect_error(tsne(exp_transform_table, group_info = test_groupInfo[, 1:3],
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  # cluster_method == 'kmedoids'
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmedoids', group_num=2, var1 = 'abc', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmedoids', group_num=2, var1 = 'euclidean', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE)
  # cluster_method == 'hclustering'
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='hclustering', group_num=2, var1 = 'abc', var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='hclustering', group_num=2, var1 = 'pearson', var2 ='abc',
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  tsne(exp_transform_table, group_info = NULL, sig_feature = NULL,
       pca=TRUE, perplexity=5, max_iter=500, cluster_method='hclustering',
       group_num=2, var1='euclidean', var2='ward.D', insert_ref_group = NULL,
       ref_group = NULL, plotly=TRUE)
  # cluster_method == 'dbscan'
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='dbscan', group_num=2, var1 = NULL, var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='dbscan', group_num=2, var1 = 1, var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  tsne(exp_transform_table, group_info = NULL, sig_feature = NULL,
       pca=TRUE, perplexity=5, max_iter=500, cluster_method='dbscan',
       group_num=2, var1=1, var2=2, insert_ref_group = NULL,
       ref_group = NULL, plotly=TRUE)
  # test Group number must be between 1 and 10
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=11, var1 = NULL, var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))
  # test Group number can not be NULL
  expect_error(tsne(exp_transform_table, group_info = NULL,
       sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
       cluster_method='kmeans', group_num=NULL, var1 = NULL, var2 = NULL,
       insert_ref_group = NULL, ref_group = NULL, plotly=TRUE))

  tsne_result <- tsne(exp_transform_table, group_info = NULL,
      sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
      cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
      insert_ref_group = NULL, ref_group = NULL, plotly=TRUE)
  expect_s3_class(tsne_result[[1]], "data.frame")
  expect_s3_class(tsne_result[[2]], "plotly")
  tsne_result <- tsne(exp_transform_table, group_info = NULL,
      sig_feature = NULL, pca=TRUE, perplexity=5, max_iter=500,
      cluster_method='kmeans', group_num=2, var1 = 'euclidean', var2 = NULL,
      insert_ref_group = NULL, ref_group = NULL, plotly=FALSE)
  expect_s3_class(tsne_result[[2]], "ggplot")
})
