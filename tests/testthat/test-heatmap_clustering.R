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
## DE - multi
data("se_multiGroup")
processed_se_multi <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_multi <- deSp_multiGroup(
    processed_se_multi, ref_group='ctrl', test='One-way ANOVA', significant='pval',
    p_cutoff=0.05, transform='log10')

deChar_se_multi <- deChar_multiGroup(
    processed_se_multi, char='class', ref_group='ctrl', post_hoc='One-way ANOVA',
    post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10')

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
    res_sp <- suppressWarnings(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='complete', type='all')
        )
    expect_heatmap(res_sp)
    res_char <- suppressWarnings(
        heatmap_clustering(deChar_se, char='Total.C', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_char)
    ## DE-multi
    res_sp_multi <- suppressWarnings(
        heatmap_clustering(deSp_se_multi, char='class', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_sp_multi)
    res_char_multi <- suppressWarnings(
        heatmap_clustering(deChar_se_multi, char='Total.C', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_char_multi)
})

test_that("heatmap_clustering handles different char parameters", {
    expect_error(
        heatmap_clustering(deSp_se, char=NULL, distfun='pearson', hclustfun='complete', type='all'),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    char_list <- list_lipid_char(processed_se)$common_list
    names(char_list) <- NULL
    for (i in seq_along(char_list)) {
        res <- suppressWarnings(
            heatmap_clustering(deSp_se, char=char_list[i], distfun='pearson', hclustfun='complete', type='all')
        )
        expect_heatmap(res)
        # multi
        res_multi <- suppressWarnings(
            heatmap_clustering(deSp_se_multi, char=char_list[i], distfun='pearson', hclustfun='complete', type='all')
        )
        expect_heatmap(res_multi)
    }
})

test_that("heatmap_clustering handles different distfun parameters", {
    for (dist in c('pearson', 'kendall', 'spearman')) {
        res <- suppressWarnings(
            heatmap_clustering(deSp_se, char='class', distfun=dist, hclustfun='complete', type='all')
        )
        expect_heatmap(res)
    }

    expect_error(
        heatmap_clustering(deSp_se, char='class', distfun='invalid', hclustfun='complete', type='all'),
        "distfun must be one of 'pearson', 'kendall', or 'spearman'."
    )
})

test_that("heatmap_clustering handles different hclustfun parameters", {
    hclust_methods <- c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
    for (method in hclust_methods) {
        res <- suppressWarnings(
            heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun=method, type='all')
        )
        expect_heatmap(res)
    }

    expect_error(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='invalid', type='all'),
        "hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'."
    )
})

test_that("heatmap_clustering handles different type parameters", {
    res_all <- suppressWarnings(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_all)

    res_sig <- suppressWarnings(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='complete', type='sig')
    )
    expect_heatmap(res_sig)

    expect_error(
        heatmap_clustering(deSp_se, char='class', distfun='pearson', hclustfun='complete', type='invalid'),
        "type must be one of 'all' or 'sig'."
    )
})

test_that("heatmap_clustering handles edge cases", {
    # Test with one row data
    feature_name <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, "feature"]
    sub <- deSp_se[SummarizedExperiment::rowData(deSp_se)$feature == feature_name, ]
    S4Vectors::metadata(sub)[["sig_deSp_result"]] <- S4Vectors::metadata(sub)[["sig_deSp_result"]][1, ]
    expect_error(
        suppressWarnings(heatmap_clustering(sub, char="class", distfun='pearson', hclustfun='complete', type='sig')),
        "Insufficient number of lipids."
    )

    # Test with no significant lipids
    no_sig_se <- deSp_se
    S4Vectors::metadata(no_sig_se)[["sig_deSp_result"]] <- "No significant lipids."
    expect_error(
        suppressWarnings(heatmap_clustering(no_sig_se, char="class", distfun='pearson', hclustfun='complete', type='sig')),
        "Insufficient number of significant lipids."
    )
})

test_that("heatmap_clustering handles different data scales", {
    # Test with all positive values
    pos_se <- deSp_se
    assay(pos_se) <- abs(assay(pos_se))
    res_pos <- suppressWarnings(
        heatmap_clustering(pos_se, char='class', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_pos)

    # Test with all negative values
    neg_se <- deSp_se
    assay(neg_se) <- -abs(assay(neg_se))
    res_neg <- suppressWarnings(
        heatmap_clustering(neg_se, char='class', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_neg)

    # Test with mixed positive and negative values
    mixed_se <- deSp_se
    res_mixed <- suppressWarnings(
        heatmap_clustering(mixed_se, char='class', distfun='pearson', hclustfun='complete', type='all')
    )
    expect_heatmap(res_mixed)
})

test_that("heatmap_clustering can handle one row data condition", {
    feature_name <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]][1, "feature"]
    sub <- deSp_se[SummarizedExperiment::rowData(deSp_se)$feature == feature_name, ]
    S4Vectors::metadata(sub)[["sig_deSp_result"]] <- S4Vectors::metadata(sub)[["sig_deSp_result"]][1, ]
    expect_error(suppressWarnings(heatmap_clustering(sub, char="class", distfun='pearson', hclustfun='complete', type='sig')),
                 "Insufficient number of lipids.")
})
