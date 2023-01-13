# Load data and create example data sets
data(profiling_exp_data)
data(profiling_lipid_char_table)

# Start tests
test_that("profiling function computed successfullly.", {
  expect_s3_class(profiling_exp_data, "data.frame")
  exp_transform_table <- data_process(profiling_exp_data,
                                      exclude_var_missing=TRUE,
                                      missing_pct_limit=50,
                                      replace_zero=TRUE, zero2what='min',
                                      xmin=0.5, replace_NA=TRUE,
                                      NA2what='min', ymin=0.5,
                                      pct_transform=TRUE,
                                      data_transform=TRUE, trans_type='log',
                                      centering=FALSE, scaling=FALSE)
  expect_s3_class(exp_transform_table, "data.frame")
  expect_equal(ncol(exp_transform_table), ncol(profiling_exp_data))

  profiling_result <- exp_profiling(profiling_exp_data)
  expect_s3_class(profiling_result$i.expr.lip, "plotly")
  expect_s3_class(profiling_result$i.p.amount, "plotly")
  expect_s3_class(profiling_result$p.hist.value, "plotly")
  ## PCA
  PCA_result <- PCA(exp_transform_table,
                    group_info = NULL, sig_feature = NULL,
                    scaling=TRUE, centering=TRUE, cluster_method='kmeans',
                    group_num=2, var1 = NULL, var2 = NULL,
                    insert_ref_group=NULL, ref_group=NULL,
                    n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
  expect_s3_class(PCA_result[[1]], "prcomp")
  expect_s3_class(PCA_result[[2]], "data.frame")
  expect_s3_class(PCA_result[[3]], "data.frame")
  expect_s3_class(PCA_result[[4]], "plotly")
  expect_s3_class(PCA_result[[5]], "plotly")
  expect_s3_class(PCA_result[[6]], "plotly")
  expect_s3_class(PCA_result[[7]], "plotly")
  ## t-SNE
  tsne_result <- tsne(exp_transform_table, group_info = NULL,
                      sig_feature = NULL, pca=TRUE, perplexity=5,
                      max_iter=500, cluster_method='kmeans',
                      group_num=2, var1 = 'euclidean', var2 = NULL,
                      insert_ref_group = NULL, ref_group = NULL, plotly=TRUE)
  expect_s3_class(tsne_result[[1]], "data.frame")
  expect_s3_class(tsne_result[[2]], "plotly")
  ## UMAP
  UMAP_result <- UMAP(exp_transform_table, group_info=NULL,
                      sig_feature=NULL, n_neighbors=15,
                      scale=TRUE, metric='euclidean', group_num=2,
                      cluster_method='kmeans', var1=NULL, var2=NULL,
                      insert_ref_group=NULL, ref_group=NULL, plotly=TRUE)
  expect_s3_class(UMAP_result[[1]], "data.frame")
  expect_s3_class(UMAP_result[[2]], "plotly")
  ## correlation heatmap
  corr_result <- corr_heatmap(exp_transform_table, corr_method="pearson",
                              distfun="maximum", hclustfun="average",
                              plotly=TRUE)
  expect_type(corr_result$sample_corr_p, "double")
  expect_type(corr_result$reorder_sample_corr_coef, "double")
  expect_type(corr_result$lipid_corr_p, "double")
  expect_type(corr_result$lipids_corr_coef, "double")
  expect_s4_class(corr_result$sample_hm, "Iheatmap")
  expect_s4_class(corr_result$lipids_hm, "Iheatmap")
  ## exp_compo_by_lipidinfo
  char_var <- colnames(profiling_lipid_char_table)[-1]
  expect_type(char_var, "character")
  compo_result <- exp_compo_by_lipidinfo(profiling_exp_data,
                                         profiling_lipid_char_table,
                                         char_var[1],
                                         plotly=TRUE)
  expect_s3_class(compo_result$p.barplot.p, "plotly")
  expect_s3_class(compo_result$p.compos, "plotly")
})
