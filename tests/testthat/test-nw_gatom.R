library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)

data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

#### Function ####
gatom_test <- function(res){
    # Verify the class of the return objects.
    expect_type(res, 'list')
    expect_s3_class(res$table_edge, 'data.frame')
    expect_s3_class(res$table_node, 'data.frame')
    expect_s3_class(res$table_reaction, 'data.frame')
    expect_s3_class(res$table_stat, 'data.frame')
    # Confirm the existence of required columns.
    expect_true(all(c('Path', 'Reaction', 'Enzyme', 'Gene') %in% colnames(res$table_reaction)))
    expect_true(all(c('Lipid', 'Log2FC', 'p-value') %in% colnames(res$table_stat)))
    expect_true(all(c(
        'from', 'to', 'from_label', 'to_label', 'label', 'dashes', 'font.size',
        'width', 'smooth', 'shadow', 'color', 'title', 'Path', 'Reaction',
        'Enzyme', 'log2FC', 'pval') %in% colnames(res$table_edge)))
    expect_true(all(c(
        'id', 'label', 'group', 'shape', 'size', 'font.size', 'shadow', 'title',
        'color', 'lipid', 'log2FC', 'pval') %in% colnames(res$table_node)))
    # Confirm whether the data frame consists of numeric values.
    expect_true(all(apply(res$table_node[,5:6], 2, function(x) all(is.numeric(x)))))
    expect_true(all(apply(res$table_edge[,c(7:8, 16:17)], 2, function(x) all(is.numeric(x)))))
    expect_true(all(apply(res$table_stat[,2:3], 2, function(x) all(is.numeric(x)))))
    # Check NA values
    expect_true(all(!is.na(res$table_reaction$Path)))
    expect_true(all(!is.na(res$table_stat$Lipid)))
    expect_true(all(apply(res$table_node[, 1:9], 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_edge[, 1:2], 2, function(x) all(!is.na(x)))))
}

#### Test ####
test_that("Test under normal conditions", {
    for(organism in c('human', 'mouse')){
        for(sp_significant in c('pval', 'padj')){
            res <- suppressWarnings(
                nw_gatom(
                    deSp_se=deSp_se_twoGroup, organism, n_lipid=50,
                    sp_significant, sp_p_cutoff=0.05, sp_FC_cutoff=1))
            gatom_test(res)
        }
    }
})

test_that("Test without significant lipids", {
    expect_warning(
        nw_gatom(
            deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
            sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1000),
        'No lipid species meet these significant cutoffs; please adjust the thresholds you have set.')
    # res <- nw_gatom(
    #     deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
    #     sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1000)
    # expect_identical(res, 'No lipid species meet these significant cutoffs; please adjust the thresholds you have set.')
})

test_that("Test with positive significant lipids", {
    S4Vectors::metadata(deSp_se_twoGroup)$all_deSp_result %<>%
        dplyr::mutate(log2FC=abs(log2FC))
    res <- suppressWarnings(
        nw_gatom(
            deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
            sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1))
    gatom_test(res)
})

test_that("Test with negative significant lipids", {
    S4Vectors::metadata(deSp_se_twoGroup)$all_deSp_result %<>%
        dplyr::mutate(log2FC=-abs(log2FC))
    res <- suppressWarnings(
        nw_gatom(
            deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
            sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1))
    gatom_test(res)
})

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
deSp_se_multi <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')

test_that("Test for incorrect input", {
    expect_error(
        suppressWarnings(
            nw_gatom(
                deSp_se=deSp_se_twoGroup, organism=NULL, n_lipid=50,
                sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1)),
        "The 'organism' parameter must be either 'human' or 'mouse'."
    )
    expect_error(
        suppressWarnings(
            nw_gatom(
                deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=0,
                sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1)),
        "The minimum value for the 'n_lipid' parameter is 1."
    )
    expect_error(
        suppressWarnings(
            nw_gatom(
                deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
                sp_significant='Pval', sp_p_cutoff=0.05, sp_FC_cutoff=1)),
        "The 'sp_significant' parameter must be either 'pval' or 'padj'."
    )
    expect_error(
        suppressWarnings(
            nw_gatom(
                deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
                sp_significant='pval', sp_p_cutoff=1.5, sp_FC_cutoff=1)),
        "The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1."
    )
    expect_error(
        suppressWarnings(
            nw_gatom(
                deSp_se=deSp_se_twoGroup, organism='mouse', n_lipid=50,
                sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=0)),
        "The minimum value for the 'sp_FC_cutoff' parameter is 1."
    )
    expect_error(suppressWarnings(
        nw_gatom(
            deSp_se=deSp_se_multi, organism='mouse', n_lipid=50,
            sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1)),
        "This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis."
    )
    sub <- processed_se_twoGroup[, SummarizedExperiment::colData(processed_se_twoGroup)$sample_name %in% c("control_01", "hfref_patient_01")]
    deSp_se <- suppressWarnings(
        deSp_twoGroup(
            sub, ref_group='ctrl', test='Wilcoxon test', significant='padj',
            p_cutoff=0.05, FC_cutoff=2, transform='log10'))
    expect_error(suppressWarnings(
        nw_gatom(
            deSp_se, organism='mouse', n_lipid=50,
            sp_significant='pval', sp_p_cutoff=0.05, sp_FC_cutoff=1)),
        "This function requires each group to have at least two samples."
    )
})
