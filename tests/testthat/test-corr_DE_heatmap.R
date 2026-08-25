library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# Load data and create example data sets
data("corr_data")
processed_se <- data_process(
    corr_data, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

data("corr_group_info")
corr_group_info[2, "Age"] <- NA
corr_data_err <- processed_se
SummarizedExperiment::colData(corr_data_err) <- S4Vectors::DataFrame(corr_group_info)

# get char list
char_list <- list_lipid_char(processed_se)$common_list
names(char_list) <- NULL

expect_corr_heatmap <- function(res) {
    if (!is.character(res)) {
        # Check if `res` is a list
        expect_true(is.list(res), info = "Output should be a list")
        expect_named(
            res,
            c("all_correlation_result", "sig_correlation_result",
              "interactive_heatmap", "static_heatmap", "heatmap_matrix"))
        expect_s3_class(res$all_correlation_result, "data.frame")
        expect_s3_class(res$sig_correlation_result, "data.frame")
        expect_s3_class(res$interactive_heatmap, "plotly")
        expect_s4_class(res$static_heatmap, "Heatmap")
        expect_true(is.matrix(res$heatmap_matrix))
    }
}

# Start tests
test_that("corr_DE_heatmap returns correct heatmap result", {
    res <- corr_DE_heatmap(
        processed_se, char=NULL, condition_col=c('Sex', 'Smoking','Exacerbations'),
        side_color_char=NULL, significant='pval',
        p_cutoff=1, adjust_p_method='BH', distfun='spearman',
        hclustfun='centroid', transform='log10', type='Sp')
    expect_corr_heatmap(res)
})

test_that("corr_DE_heatmap can handle error inputs", {
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char='char', condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Char'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='p', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=1.2,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='other', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='Spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='Centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log2', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='other'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char="class", condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char='char', significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
            side_color_char='char', significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Char'), error=TRUE)
})

test_that("corr_DE_heatmap can handle condition_col and adjusted_col inputs.", {
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c("ABC", "Emphysema"),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=c("Emphysema"),
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
    expect_snapshot(
        corr_DE_heatmap(
            processed_se, char=NULL, condition_col=NULL,
            side_color_char=NULL, significant='pval', p_cutoff=0.05,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp'), error=TRUE)
})

test_that("corr_DE_heatmap can handle all adjust_p_method.", {
    adjust_p_list <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')
    for (adjust_p in adjust_p_list) {
        res_padj <- corr_DE_heatmap(
            processed_se, char=NULL,
            condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='padj', p_cutoff=1,
            adjust_p_method=adjust_p, distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp')
        expect_corr_heatmap(res_padj)
    }
})

test_that("corr_DE_heatmap can handle all distfun & hclustfun.", {
    distfun_list <- c('pearson', 'spearman', 'kendall')
    hclustfun_list <- c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
    for (distfun in distfun_list) {
        for (hclustfun in hclustfun_list) {
            res <- corr_DE_heatmap(
                processed_se, char=NULL,
                condition_col=c('Sex', 'Smoking'),
                side_color_char=NULL, significant='pval', p_cutoff=1,
                adjust_p_method='BH', distfun=distfun, hclustfun=hclustfun,
                transform='log10', type='Sp')
            expect_corr_heatmap(res)
        }
    }
})

test_that("corr_DE_heatmap can handle all transform.", {
    transform_list <- c('none', 'log10', 'square', 'cube')
    for (transform in transform_list) {
        res <- corr_DE_heatmap(
            processed_se, char=NULL,
            condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=1,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform=transform, type='Sp')
        expect_corr_heatmap(res)
    }
})

test_that("corr_DE_heatmap can handle all lipid characteristics. (Sp)", {
    for (char in char_list) {
        res <- corr_DE_heatmap(
            processed_se, char=NULL,
            condition_col=c('Sex', 'Smoking'),
            side_color_char=char, significant='pval', p_cutoff=1,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Sp')
        expect_corr_heatmap(res)
    }
})

test_that("corr_DE_heatmap can handle all lipid characteristics. (Char)", {
    char_seq <- char_list[
        !(char_list %in% c("Category", "FA.OH", "FA.Unsaturation.Category1", "Total.OH", "Intrinsic.Curvature"))]
    for (char in char_seq) {
        res <- corr_DE_heatmap(
            processed_se, char=char,
            condition_col=c('Sex', 'Smoking'),
            side_color_char=NULL, significant='pval', p_cutoff=1,
            adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
            transform='log10', type='Char')
        expect_corr_heatmap(res)
    }
})
