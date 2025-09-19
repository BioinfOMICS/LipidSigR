library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# load example data
data("profiling_data")
profiling_se <- data_process(
    profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

## DE
data("de_data_twoGroup")
processed_se <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se <- deSp_twoGroup(
    processed_se, ref_group='ctrl', test='t-test', significant='padj',
    p_cutoff=0.05, FC_cutoff=2, transform='log10')
deChar_se <- deChar_twoGroup(
    processed_se, char="Total.FA", ref_group="ctrl", test='t-test',
    significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')

## DE - multi
data("se_multiGroup")
processed_se_multi <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_multi <- deSp_multiGroup(
    processed_se_multi, ref_group='ctrl', test='One-way ANOVA', significant='pval',
    p_cutoff=0.05, transform='log10')

deChar_se_multi <- deChar_multiGroup(
    processed_se_multi, char='class', ref_group='ctrl', post_hoc='One-way ANOVA',
    post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10')


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

### PCA ---------------------
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
    ## DE-multi
    PCA <- suppressWarnings(
        dr_pca(deSp_se_multi, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=2)
    )
    expect_PCA(PCA)
    PCA <- suppressWarnings(
        dr_pca(deChar_se_multi, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=2)
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
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=NULL, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "scaling must be a logical value.")
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering="yes",
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "centering must be a logical value.")
})
# Test for clustering
test_that("dr_pca clustering & related parameters", {
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='k-means', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "clustering must be one of 'kmeans', 'kmedoids', 'hclustering', 'dbscan', or 'group_info'.")
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=0, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "cluster_num must be a numeric value between ")
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=10, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10))
    expect_PCA(PCA)
    # kmedoids
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmedoids', cluster_num=2, kmedoids_metric="euclidean", distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
    )
    expect_PCA(PCA)
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmedoids', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        'kmedoids_metric must be one of "euclidean" or "manhattan".')
    # hclustering
    for (distfun in c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
        for (hclustfun in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
            PCA <- suppressWarnings(
                dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
                       clustering='hclustering', cluster_num=2, kmedoids_metric=NULL, distfun,
                       hclustfun, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)
            )
            expect_PCA(PCA)
        }
    }
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='hclustering', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun="ward.D", eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        'distfun must be one of "pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski".')
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='hclustering', cluster_num=2, kmedoids_metric=NULL, distfun="maximum",
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        'hclustfun must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid".')
    # dbscan
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='dbscan', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=10, minPts=100, feature_contrib_pc=c(1,2), plot_topN=10)
    )
    expect_PCA(PCA)
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='dbscan', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=0, minPts=0, feature_contrib_pc=c(1,2), plot_topN=10)
    )
    expect_PCA(PCA)
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='dbscan', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=100, feature_contrib_pc=c(1,2), plot_topN=10)),
        'eps must be must be numeric value >= 0.')
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='dbscan', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=0, minPts=-1, feature_contrib_pc=c(1,2), plot_topN=10)),
        'minPts must be must be numeric value >= 0.')
})
# Test plot_topN
test_that("dr_pca plot_topN parameters", {
    PCA <- suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=9.5)
    )
    expect_PCA(PCA)
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=0)),
        "plot_topN must be a numeric value less than the number of lipids.")
})
# Test for Insufficient number of significant lipids
test_that("dr_pca handles insufficient number of significant lipids condition", {
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, ]
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "Insufficient number of significant lipids.")
    S4Vectors::metadata(deSp_se)[["sig_deSp_result"]] <- "No significant lipid."
    expect_error(suppressWarnings(
        dr_pca(deSp_se, scaling=TRUE, centering=TRUE,
               clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
               hclustfun=NULL, eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)),
        "Insufficient number of significant lipids.")
})


### tsne -----------------
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
    ## DE-multi
    tsne <- suppressWarnings(
        dr_tsne(deSp_se_multi, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
    tsne <- suppressWarnings(
        dr_tsne(deChar_se_multi, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
})
# Test for error parameter
test_that("dr_tsne can return error message for incorrect input", {
    expect_error(suppressWarnings(
        dr_tsne(profiling_se, pca=NULL, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)),
        "pca must be a logical value.")
    tsne <- suppressWarnings(
        dr_tsne(deSp_se, pca=FALSE, perplexity=5, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
    expect_error(suppressWarnings(
        dr_tsne(profiling_se, pca=TRUE, perplexity=-1, max_iter=500, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)),
        "perplexity cannot be a negative value.")
    expect_error(suppressWarnings(
        dr_tsne(profiling_se, pca=TRUE, perplexity=5, max_iter=-1, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)),
        "max_iter cannot be a negative value.")
    tsne <- suppressWarnings(
        dr_tsne(deSp_se, pca=TRUE, perplexity=5, max_iter=1, clustering='kmeans',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_tsne(tsne)
})
# Test for clustering
test_that("dr_tsne dbscan clustering", {
    tsne <- suppressWarnings(
        dr_tsne(deSp_se, pca=TRUE, perplexity=5, max_iter=500, clustering='dbscan',
                cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=0, minPts=1)
    )
    expect_tsne(tsne)
})

# UMAP-----------------------
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
    ## DE-multi
    umap <- suppressWarnings(
        dr_umap(deSp_se_multi, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deChar_se_multi, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
})
# Test for incorrect input
test_that("dr_umap returns error message for incorrect input", {
    expect_error(
        suppressWarnings(
            dr_umap(profiling_se, n_neighbors=0, scaling=TRUE, umap_metric='euclidean',
                    clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                    hclustfun=NULL, eps=NULL, minPts=NULL)),
        "n_neighbors must be >= 2.")
    umap <- suppressWarnings(
        dr_umap(profiling_se, n_neighbors=2, scaling=TRUE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL))
    expect_umap(umap)
    # scaling
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="none", umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling=FALSE, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling=NULL, umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="Z", umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="scale", umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="range", umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="colrange", umap_metric='euclidean',
                clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=NULL, minPts=NULL)
    )
    expect_umap(umap)
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling="maxabs", umap_metric='euclidean',
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

    expect_error(
        suppressWarnings(
            dr_umap(profiling_se, n_neighbors=150, scaling="YES", umap_metric='euclidean',
                    clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                    hclustfun=NULL, eps=NULL, minPts=NULL)),
        'scaling must be one of "none", FALSE, NULL, "Z", "scale", TRUE, "maxabs", "range", or "colrange".')
    # umap_metric
    for (umap_metric in c('euclidean', 'cosine', 'manhattan', 'hamming', "correlation")) {
        umap <- suppressWarnings(
            dr_umap(deChar_se_multi, n_neighbors=15, scaling=TRUE, umap_metric,
                    clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                    hclustfun=NULL, eps=NULL, minPts=NULL)
        )
        expect_umap(umap)
    }
    expect_error(
        suppressWarnings(
            dr_umap(profiling_se, n_neighbors=150, scaling=TRUE, umap_metric='Euclidean',
                    clustering='kmeans', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                    hclustfun=NULL, eps=NULL, minPts=NULL)),
        "umap_metric must be one of 'euclidean', 'cosine', 'manhattan', 'hamming', or 'correlation'.")
})
# Test for clustering
test_that("dr_umap dbscan clustering", {
    umap <- suppressWarnings(
        dr_umap(deSp_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
                clustering='dbscan', cluster_num=2, kmedoids_metric=NULL, distfun=NULL,
                hclustfun=NULL, eps=0, minPts=10)
    )
    expect_umap(umap)
})


# PLSDA ---------------------
# Test for PLSDA function and output structure
test_that("dr_plsda returns correct output results", {
    expect_error(
        dr_plsda(
            profiling_se, ncomp=2, scaling=TRUE, clustering='group_info',
            cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL),
        "Incorrect SummarizedExperiment structure.")
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
    ## DE-multi
    plsda <- suppressWarnings(
        dr_plsda(deSp_se_multi, ncomp=2, scaling=TRUE, clustering='group_info',
                 cluster_num=3, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)
    )
    expect_plsda(plsda)
    plsda <- suppressWarnings(
        dr_plsda(deChar_se_multi, ncomp=2, scaling=TRUE, clustering='group_info',
                 cluster_num=3, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)
    )
    expect_plsda(plsda)
})
# Test incorrect parameter
test_that("dr_plsda returns error message for incorrect input.", {
    expect_error(suppressWarnings(
        dr_plsda(deSp_se, ncomp=-1, scaling=TRUE, clustering='group_info',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)),
        "ncomp cannot be a negative value or NULL.")
    expect_error(suppressWarnings(
        dr_plsda(deSp_se, ncomp=2, scaling=1, clustering='group_info',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)),
        "scaling must be a logical value.")
    plsda <- suppressWarnings(
        dr_plsda(deChar_se_multi, ncomp=2, scaling=FALSE, clustering='group_info',
                 cluster_num=3, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=NULL, minPts=NULL)
    )
    expect_plsda(plsda)
})
# Test for clustering
test_that("dr_plsda dbscan clustering", {
    plsda <- suppressWarnings(
        dr_plsda(deChar_se_multi, ncomp=2, scaling=FALSE, clustering='dbscan',
                 cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
                 eps=0, minPts=0)
    )
    expect_plsda(plsda)
})
