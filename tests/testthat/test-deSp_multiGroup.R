library(testthat)
#library(SummarizedExperiment)
data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
p_cutoff <- 0.05

#### Function ####
DE_multi_se_format <- function(se) {
    # Check if `se` is an SE object
    expect_true(inherits(se, "SummarizedExperiment"),
                info = "Object is not an SummarizedExperiment")
    # Check if mandatory components are present
    expect_true("assays" %in% slotNames(se),
                info = "Assays slot is missing")
    # expect_true("rowData" %in% slotNames(se),
    #             info = "rowData slot is missing")
    expect_true("colData" %in% slotNames(se),
                info = "colData slot is missing")
    expect_true("metadata" %in% slotNames(se),
                info = "metadata slot is missing")
    # Check the structure of assays, rowData, and colData
    # Check expression data exists and is a matrix
    expect_true(is.matrix(SummarizedExperiment::assay(se)),
                info = "Expression data is not a matrix")
    expect_true(all(dim(SummarizedExperiment::assay(se)) > 0),
                info = "Expression data has non-positive dimensions")
}
DE_multi_sig_lipid <- function(DE.res, p_cutoff){
    # Verify if the significant lipids indeed represent significance.
    expect_equal(DE.res$feature[which(DE.res$pval < p_cutoff)],
                 DE.res$feature[which(DE.res$sig_pval == 'yes')])
    expect_equal(DE.res$feature[which(DE.res$padj < p_cutoff)],
                 DE.res$feature[which(DE.res$sig_padj == 'yes')])
}
DE_multi_result_table <- function(DE.res){
    # Verify if it is a data frame.
    expect_true(is.data.frame(DE.res))
    # Confirm whether columns that cannot be NA indeed do not contain any NA values.
    expect_true(all(!is.na(DE.res$feature)))
    # Confirm the existence of required columns.
    expect_true(all(c('feature', 'pval', 'negLog10pval', 'padj', 'negLog10padj',
                      'sig_pval', 'sig_padj') %in% colnames(DE.res)))
    expect_false('log2FC' %in% colnames(DE.res))
}
DE_multi_test <- function(deSp_se, significant, p_cutoff){
    # Verify the SE object returned.
    DE_multi_se_format(se=deSp_se)

    abundance <- .extract_df(deSp_se, type="abundance")
    all_deSp_result <- S4Vectors::metadata(deSp_se)$all_deSp_result
    sig_deSp_result <- S4Vectors::metadata(deSp_se)$sig_deSp_result

    # Confirm if the lipid features are the same.
    expect_equal(sort(abundance$feature), sort(all_deSp_result$feature))
    # Confirm the correctness of significance.
    DE_multi_sig_lipid(DE.res=all_deSp_result, p_cutoff=p_cutoff)
    DE_multi_sig_lipid(DE.res=sig_deSp_result, p_cutoff=p_cutoff)
    # Confirm the existence of required columns.
    DE_multi_result_table(DE.res=all_deSp_result)
    DE_multi_result_table(DE.res=sig_deSp_result)
    # Confirm the correctness of significance and the top 20 lipids
    if(significant == 'pval'){
        expect_true(all(sig_deSp_result$sig_pval == 'yes'))
    }else if(significant == 'padj'){
        expect_true(all(sig_deSp_result$sig_padj == 'yes'))
    }
    # Verify if it is a data frame.
    expect_true(is.data.frame(S4Vectors::metadata(deSp_se)$processed_abundance))
    expect_true(is.character(S4Vectors::metadata(deSp_se)$significant))
    expect_true(is.numeric(S4Vectors::metadata(deSp_se)$p_cutoff))
    expect_true(is.character(S4Vectors::metadata(deSp_se)$transform))
}

#### Test ####
test_that("Test under normal conditions", {
    for(test in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                deSp_se <- deSp_multiGroup(
                    processed_se=processed_se_multiGroup, ref_group='ctrl', test,
                    significant, p_cutoff, transform)
                DE_multi_test(deSp_se, significant, p_cutoff)
            }
        }
    }
})






