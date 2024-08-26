library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# load example data
data("de_data_twoGroup")
processed_se <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage')
deSp_se <- deSp_twoGroup(
    processed_se, ref_group='ctrl', test='t-test',
    significant='padj', p_cutoff=0.05, FC_cutoff=2, transform='log10')
deChar_se <- deChar_twoGroup(
    processed_se, char="Total.C", ref_group="ctrl", test='t-test',
    significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')

expect_heatmap <- function(res) {
    # Check if `res` is a list
    expect_true(is.list(res), info = "Output should be a list")
    expect_true("interactive_heatmap" %in% names(res), info = "'interactive_heatmap' should be in the list")
    expect_true("static_heatmap" %in% names(res), info = "'static_heatmap' should be in the list")
    expect_true("corr_coef_matrix" %in% names(res), info = "'corr_coef_matrix' should be in the list")
    expect_s4_class(res$interactive_heatmap, "IheatmapHorizontal")
    expect_true(is.matrix(res$corr_coef_matrix))
    expect_s3_class(res$static_heatmap, "recordedplot")
}

test_that("heatmap_clustering function outputs correct structure", {
    res <- suppressWarnings(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='complete', type='all')
        )
    expect_heatmap(res)
    # res <- suppressWarnings(
    #     heatmap_clustering(deChar_se, char='Total.C', distfun='pearson', hclustfun='complete', type='all')
    # )
    # expect_heatmap(res)
})

test_that("heatmap_clustering char parameter", {
    expect_error(
        heatmap_clustering(deSp_se, char=NULL, distfun='pearson', hclustfun='complete', type='all'),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    # char_list <- list_lipid_char(processed_se)$common_list
    # names(char_list) <- NULL
    # for (i in seq(char_list)) {
    #     res <- suppressWarnings(
    #         heatmap_clustering(deSp_se, char=char_list[i], distfun='pearson', hclustfun='complete', type='all')
    #     )
    #     expect_heatmap(res)
    # }
})

test_that("heatmap_clustering can handle one row data condition", {
    feature_name <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, "feature"]
    sub <- deSp_se[SummarizedExperiment::rowData(deSp_se)$feature == feature_name, ]
    S4Vectors::metadata(sub)[["sig_deSp_result"]] <- S4Vectors::metadata(sub)[["sig_deSp_result"]][1, ]
    expect_error(suppressWarnings(heatmap_clustering(sub, char="class", distfun='pearson', hclustfun='complete', type='sig')),
                 "Insufficient number of lipids.")
})
