library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)
library(tidyr)

data("se_multiGroup")
processed_se <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
post_hoc_p_cutoff <- 0.05

data("de_data_twoGroup")
processed_se_two <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')

#### Function ####
DE_multi_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info="Object is not an SummarizedExperiment")
    # Check if mandatory components are present
    expect_true("assays" %in% slotNames(se),
                info="Assays slot is missing")
    # expect_true("rowData" %in% slotNames(se),
    #             info="rowData slot is missing")
    expect_true("colData" %in% slotNames(se),
                info="colData slot is missing")
    expect_true("metadata" %in% slotNames(se),
                info="metadata slot is missing")
    # Check the structure of assays, rowData, and colData
    # Check expression data exists and is a matrix
    expect_true(is.matrix(SummarizedExperiment::assay(se)),
                info="Expression data is not a matrix")
    expect_true(all(dim(SummarizedExperiment::assay(se)) > 0),
                info="Expression data has non-positive dimensions")
}
DE_multi_sig_lipid <- function(DE.res, post_hoc_p_cutoff){
    # Verify if the significant lipids indeed represent significance.
    expect_equal(DE.res$feature[which(DE.res$post_hoc_pval < post_hoc_p_cutoff)],
                 DE.res$feature[which(DE.res$post_hoc_sig_pval == 'yes')])
    expect_equal(DE.res$feature[which(DE.res$post_hoc_padj < post_hoc_p_cutoff)],
                 DE.res$feature[which(DE.res$post_hoc_sig_padj == 'yes')])
}
DE_multi_result_table <- function(DE.res){
    # Verify if it is a data frame.
    expect_true(is.data.frame(DE.res))
    # Confirm whether columns that cannot be NA indeed do not contain any NA values.
    expect_true(all(!is.na(DE.res$feature)))
    # Confirm the existence of required columns.
    expect_true(all(c(
        'characteristic', 'feature', 'post_hoc_pval', 'post_hoc_negLog10pval',
        'post_hoc_padj', 'post_hoc_negLog10padj', 'post_hoc_sig_pval',
        'post_hoc_sig_padj') %in% colnames(DE.res)))
    expect_false('log2FC' %in% colnames(DE.res))
}
DE_char_multi_test <- function(subChar_se, post_hoc_sig, post_hoc_p_cutoff){
    # Verify the SE object returned.
    DE_multi_se_format(se=subChar_se)

    abundance <- .extract_df(subChar_se, type="abundance")
    all_deChar_result <- S4Vectors::metadata(subChar_se)$all_deChar_result
    sig_deChar_result <- S4Vectors::metadata(subChar_se)$sig_deChar_result

    # Confirm if the lipid features are the same.
    expect_equal(sort(abundance$feature),
                 sort(paste0(all_deChar_result$sub_feature, '_', all_deChar_result$feature)))
    # Confirm the correctness of significance.
    DE_multi_sig_lipid(DE.res=all_deChar_result, post_hoc_p_cutoff=post_hoc_p_cutoff)
    DE_multi_sig_lipid(DE.res=sig_deChar_result, post_hoc_p_cutoff=post_hoc_p_cutoff)
    # Confirm the existence of required columns.
    DE_multi_result_table(DE.res=all_deChar_result)
    DE_multi_result_table(DE.res=sig_deChar_result)
    # Confirm the correctness of significance and the top 20 lipids
    if(post_hoc_sig == 'pval'){
        expect_true(all(sig_deChar_result$post_hoc_sig_pval == 'yes'))
    }else if(post_hoc_sig == 'padj'){
        expect_true(all(sig_deChar_result$post_hoc_sig_padj == 'yes'))
    }
    # Verify if it is a data frame.
    expect_true(is.data.frame(S4Vectors::metadata(subChar_se)$processed_abundance))
    expect_true(is.character(S4Vectors::metadata(subChar_se)$post_hoc_sig))
    expect_true(is.numeric(S4Vectors::metadata(subChar_se)$post_hoc_p_cutoff))
}

#### Test ####
test_that("Test under normal conditions: Lipid characteristics", {
    ## General lipid characteristics
    for(post_hoc in c('One-way ANOVA', 'Kruskal-Wallis test')){
        for(post_hoc_sig in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                subChar_se <- subChar_multiGroup(
                    processed_se, char='Total.C', subChar='class', ref_group='ctrl',
                    post_hoc, post_hoc_sig, post_hoc_p_cutoff, transform)
                DE_char_multi_test(subChar_se, post_hoc_sig, post_hoc_p_cutoff)
            }
        }
    }
})

# Test subChar_multiGroup function
test_that("subChar_multiGroup can handle invalid inputs.", {
    # Test error handling
    expect_error(
        subChar_multiGroup(
            processed_se, char='Invalid', subChar='class', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "The 'char' parameter must be a numeric lipid characteristic, such as 'Total.C'. Please execute list_lipid_char function to view the valid lipid characteristics you can input."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='Total.C', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "The 'char' and 'subChar' parameters cannot be the same."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='Invalid', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "The 'subChar' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='class', ref_group='Invalid',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "ref_group must be one of the group names in the 'group' column of the group information table."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='class', ref_group='ctrl', post_hoc='Invalid',
            post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "The 'post_hoc' parameter must be either 'One-way ANOVA' or 'Kruskal-Wallis test'."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='class', ref_group='ctrl', post_hoc='One-way ANOVA',
            post_hoc_sig='Pvalue', post_hoc_p_cutoff=0.05, transform='log10'),
        "The 'post_hoc_sig' parameter must be either 'pval' or 'padj'."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='class', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=2, transform='log10'),
        "The 'post_hoc_p_cutoff' parameter must be a numeric value between 0 and 1."
        )
    expect_error(
        subChar_multiGroup(
            processed_se, char='Total.C', subChar='class', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log2'),
        "The 'transform' parameter must be one of 'none', 'log10','square', or 'cube'."
        )
    expect_error(
        subChar_multiGroup(
            processed_se_two, char='Total.C', subChar='class', ref_group='ctrl',
            post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10'),
        "This function is used for multiple groups. If your data consists of two groups, please use subChar_twoGroup function for analysis."
    )
})

# Test helper functions
test_that("Helper functions work correctly", {
    subchar_result <- subChar_multiGroup(
        processed_se=processed_se, char='Total.C', subChar='class',
        ref_group='ctrl', post_hoc='One-way ANOVA', post_hoc_sig='pval',
        post_hoc_p_cutoff=0.05, transform='log10'
    )

    # Test DE_multi_se_format
    expect_no_error(DE_multi_se_format(subchar_result))

    # Test DE_multi_sig_lipid
    all_result <- S4Vectors::metadata(subchar_result)$all_deChar_result
    expect_no_error(DE_multi_sig_lipid(all_result, 0.05))

    # Test DE_multi_result_table
    expect_no_error(DE_multi_result_table(all_result))

    # Test DE_char_multi_test
    expect_no_error(DE_char_multi_test(subchar_result, 'pval', 0.05))
})

# Test edge cases
test_that("Edge cases are handled correctly", {
    # Test with a very high p-value cutoff
    result_high_p <- subChar_multiGroup(
        processed_se=processed_se, char='Total.C', subChar='class',
        ref_group='ctrl', post_hoc='One-way ANOVA', post_hoc_sig='pval',
        post_hoc_p_cutoff=0.99, transform='log10'
    )
    DE_char_multi_test(subChar_se=result_high_p, post_hoc_sig='pval', post_hoc_p_cutoff=0.99)

    # Test with a very low p-value cutoff
    result_low_p <- subChar_multiGroup(
        processed_se=processed_se, char='Total.C', subChar='class',
        ref_group='ctrl', post_hoc='One-way ANOVA', post_hoc_sig='pval',
        post_hoc_p_cutoff=1e-10, transform='log10'
    )
    DE_char_multi_test(subChar_se=result_low_p, post_hoc_sig='pval', post_hoc_p_cutoff=1e-10)
})

