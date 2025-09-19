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
corr_group_info[2, 3] <- NA
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
        expect_s3_class(res$static_heatmap, "recordedplot")
        expect_true(is.matrix(res$heatmap_matrix))
    }
}

# Start tests
test_that("corr_cor_heatmap returns correct heatmap result", {
    res <- corr_cor_heatmap(
        processed_se, char=NULL,
        condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
        side_color_char=NULL, correlation='pearson', significant='pval',
        p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
        distfun='spearman', hclustfun='average', heatmap_col='statistic',
        transform='log10', type='Sp')
    expect_corr_heatmap(res)
})

test_that("corr_cor_heatmap can handle all error inputs.", {
    expect_error(
        corr_cor_heatmap(
            processed_se, char="char",
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Char'),
        "Wrong char input, you can view the available char list by list_lipid_char function."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation=NULL, significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "correlation must be one of 'pearson', 'spearman', or 'kendall'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='p',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "significant must be one of 'pval' or 'padj'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=NULL, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "p_cutoff must be a numeric value between 0 and 1."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method=NULL, cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "adjust_p_method must be one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=2,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "cor_coef_cutoff must be a numeric value between 0 and 1."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun=NULL, hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "distfun must be one of 'pearson', 'spearman', or 'kendall'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='Average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='t_statistic',
            transform='log10', type='Sp'),
        "heatmap_col must be 'cor_coef' or 'statistic'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log2', type='Sp'),
        "transform must be one of 'none', 'log10', 'square', or 'cube'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type=NULL),
        "type must be one of 'Sp' or 'Char'."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char='class',
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char='class', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Char'),
        "If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char='char', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function."
    )
})

test_that("corr_cor_heatmap can handle all adjust_p_method.", {
    adjust_p_list <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')
    for (adjust_p in adjust_p_list) {
        res_padj <- corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method=adjust_p, cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp')
        expect_corr_heatmap(res_padj)
    }
})

test_that("corr_cor_heatmap can handle condition_col and adjusted_col inputs.", {
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("ABC", "Emphysema", "Exacerbations"),
            side_color_char='class', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "The condition_col must be column names selected from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("Exacerbations"),
            side_color_char='class', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "The condition_col must include at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1."
    )
    expect_error(
        corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=NULL,
            side_color_char='class', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "The condition_col parameter cannot be null. Select at least two columns from FEV1_FVC, Emphysema, Exacerbations, Age, Sex, Smoking, BMI, FEV1."
    )
    expect_error(
        corr_cor_heatmap(
            corr_data_err, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char='class', correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp'),
        "Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values."
    )
})

test_that("corr_cor_heatmap can handle all distfun & hclustfun.", {
    distfun_list <- c('pearson', 'spearman', 'kendall')
    hclustfun_list <- c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
    for (distfun in distfun_list) {
        for (hclustfun in hclustfun_list) {
            res <- corr_cor_heatmap(
                processed_se, char=NULL,
                condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
                side_color_char=NULL, correlation='pearson', significant='pval',
                p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
                distfun=distfun, hclustfun=hclustfun, heatmap_col='statistic',
                transform='log10', type='Sp')
            expect_corr_heatmap(res)
        }
    }
})

test_that("corr_cor_heatmap can handle all heatmap_col & transform.", {
    transform_list <- c('none', 'log10', 'square', 'cube')
    for (transform in transform_list) {
        res_cor <- corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='cor_coef',
            transform=transform, type='Sp')
        expect_corr_heatmap(res_cor)
        res <- corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform=transform, type='Sp')
        expect_corr_heatmap(res)
    }
})

test_that("corr_cor_heatmap can handle all lipid characteristics (Sp).", {
    for (char in char_list) {
        res <- corr_cor_heatmap(
            processed_se, char=NULL,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=char, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Sp')
        expect_corr_heatmap(res)
    }
})

test_that("corr_cor_heatmap can handle all lipid characteristics (Char).", {
    char_seq <- char_list[
        !(char_list %in% c("Category", "FA.OH", "FA.Unsaturation.Category1", "Total.OH", "Intrinsic.Curvature"))]
    for (char in char_seq) {
        res <- corr_cor_heatmap(
            processed_se, char=char,
            condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
            side_color_char=NULL, correlation='pearson', significant='pval',
            p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0,
            distfun='spearman', hclustfun='average', heatmap_col='statistic',
            transform='log10', type='Char')
        expect_corr_heatmap(res)
    }
})
