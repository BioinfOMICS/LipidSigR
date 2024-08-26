library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# load example data
data("profiling_data")
profiling_se <- data_process(
    profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

## DE
data("de_data_twoGroup")
processed_se <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se <- deSp_twoGroup(
    processed_se, ref_group='ctrl', test='t-test', significant='padj',
    p_cutoff=0.05, FC_cutoff=2, transform='log10')
deChar_se <- deChar_twoGroup(
    processed_se, char="Total.FA", ref_group="ctrl", test='t-test',
    significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')

expect_PCA <- function(res){
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("pca_rotated_data" %in% names(res), info = "'pca_rotated_data' should be in the list")
    expect_true("table_pca_contribution" %in% names(res), info = "'table_pca_contribution' should be in the list")
    expect_true("interactive_pca" %in% names(res), info = "'interactive_pca' should be in the list")
    expect_true("static_pca" %in% names(res), info = "'static_pca' should be in the list")
    expect_true("interactive_screePlot" %in% names(res), info = "'interactive_screePlot' should be in the list")
    expect_true("static_screePlot" %in% names(res), info = "'static_screePlot' should be in the list")
    expect_true("interactive_feature_contribution" %in% names(res), info = "'interactive_feature_contribution' should be in the list")
    expect_true("static_feature_contribution" %in% names(res), info = "'static_feature_contribution' should be in the list")
    expect_true("interactive_variablePlot" %in% names(res), info = "'interactive_variablePlot' should be in the list")
    expect_true("static_variablePlot" %in% names(res), info = "'static_variablePlot' should be in the list")

    expect_s3_class(res$pca_rotated_data, "data.frame")
    expect_s3_class(res$table_pca_contribution, "data.frame")
    # Check print_static parameter
    expect_s3_class(res$interactive_pca, "plotly")
    expect_s3_class(res$interactive_screePlot, "plotly")
    expect_s3_class(res$interactive_feature_contribution, "plotly")
    expect_s3_class(res$interactive_variablePlot, "plotly")
    expect_s3_class(res$static_pca, "ggplot")
    expect_s3_class(res$static_screePlot, "ggplot")
    expect_s3_class(res$static_feature_contribution, "ggplot")
    expect_s3_class(res$static_variablePlot, "ggplot")
}

expect_tsne <- function(res){
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("tsne_result" %in% names(res), info = "'tsne_result' should be in the list")
    expect_true("interactive_tsne" %in% names(res), info = "'interactive_tsne' should be in the list")
    expect_true("static_tsne" %in% names(res), info = "'static_tsne' should be in the list")
    expect_s3_class(res$tsne_result, "data.frame")
    # Check print_static parameter
    expect_s3_class(res$interactive_tsne, "plotly")
    expect_s3_class(res$static_tsne, "ggplot")
}

expect_umap <- function(res){
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("umap_result" %in% names(res), info = "'umap_result' should be in the list")
    expect_true("interactive_umap" %in% names(res), info = "'interactive_umap' should be in the list")
    expect_true("static_umap" %in% names(res), info = "'static_umap' should be in the list")
    expect_s3_class(res$umap_result, "data.frame")
    # Check print_static parameter
    expect_s3_class(res$interactive_umap, "plotly")
    expect_s3_class(res$static_umap, "ggplot")
}

expect_plsda <- function(res){
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("plsda_result" %in% names(res), info = "'plsda_result' should be in the list")
    expect_true("table_plsda_loading" %in% names(res), info = "'table_plsda_loading' should be in the list")
    expect_true("interacitve_plsda" %in% names(res), info = "'interacitve_plsda' should be in the list")
    expect_true("static_plsda" %in% names(res), info = "'static_plsda' should be in the list")
    expect_true("interactive_loadingPlot" %in% names(res), info = "'interactive_loadingPlot' should be in the list")
    expect_true("static_loadingPlot" %in% names(res), info = "'static_loadingPlot' should be in the list")
    expect_s3_class(res$plsda_result, "data.frame")
    expect_s3_class(res$table_plsda_loading, "data.frame")
    # Check print_static parameter
    expect_s3_class(res$interacitve_plsda, "plotly")
    expect_s3_class(res$interactive_loadingPlot, "plotly")
    expect_s3_class(res$static_plsda, "ggplot")
    expect_s3_class(res$static_loadingPlot, "ggplot")
}


# Test for PCA function and output structure
test_that("dr_pca returns correct output results", {
    PCA <- suppressWarnings(
        dr_pca(profiling_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
    )
    expect_PCA(PCA)
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
    )
    expect_PCA(PCA)
    PCA <- suppressWarnings(
        dr_pca(deChar_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
        )
    expect_PCA(PCA)
})

# Test scaling & centering in PCA
test_that("dr_pca scaling & centering parameters", {
    PCA <- suppressWarnings(
        dr_pca(deChar_se, scaling=FALSE, centering=FALSE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
        )
    expect_PCA(PCA)
})

# Test for tsne function and output structure
test_that("dr_tsne returns correct output results", {
    tsne <- suppressWarnings(
        dr_tsne(profiling_se, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
        )
    expect_tsne(tsne)
    tsne <- suppressWarnings(
        dr_tsne(deSp_se, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
    tsne <- suppressWarnings(
        dr_tsne(deChar_se, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
})

# Test for UMAP function and output structure
test_that("dr_umap returns correct output results", {
    umap <- suppressWarnings(
        dr_umap(profiling_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deChar_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
})

# Test for PLSDA function and output structure
test_that("dr_plsda returns correct output results", {
    expect_error(
        dr_plsda(profiling_se, ncomp=2, scaling=TRUE, clustering='group_info',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL), "Insufficient number of significant lipids.")
    plsda <- suppressWarnings(
        dr_plsda(deSp_se, ncomp=2, scaling=TRUE, clustering='group_info',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)
    )
    expect_plsda(plsda)
    plsda <- suppressWarnings(
        dr_plsda(deChar_se, ncomp=2, scaling=TRUE, clustering='group_info',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)
    )
    expect_plsda(plsda)
})
