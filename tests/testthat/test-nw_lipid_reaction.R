library(testthat)
library(SummarizedExperiment)
library(tidyverse)

data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

#### Function ####
lipid_reaction_test <- function(res){
    # Verify the class of the return objects.
    expect_type(res, 'list')
    expect_s3_class(res$table_edge, 'data.frame')
    expect_s3_class(res$table_node, 'data.frame')
    expect_s3_class(res$table_reaction, 'data.frame')
    expect_s3_class(res$table_stat, 'data.frame')
    # Confirm the existence of required columns.
    expect_true(all(c('Path', 'Reaction', 'Enzyme', 'Gene') %in% colnames(res$table_reaction)))
    expect_true(all(c('Level', 'Lipid', 'Log2FC', 'Significance') %in% colnames(res$table_stat)))
    expect_true(all(c('from', 'to', 'label', 'title', 'color', 'dashes', 'font.size',
                      'width', 'arrows', 'smooth', 'shadow')
                    %in% colnames(res$table_edge)))
    expect_true(all(c('id', 'label', 'group', 'shape', 'size', 'font.size',
                      'shadow', 'title', 'color') %in% colnames(res$table_node)))
    # Confirm whether the data frame consists of numeric values.
    expect_true(all(apply(res$table_node[,5:6], 2, function(x) all(is.numeric(x)))))
    expect_true(all(apply(res$table_edge[,7:8], 2, function(x) all(is.numeric(x)))))
    expect_true(is.numeric(res$table_stat$Log2FC))
    # Check NA values
    expect_true(all(!is.na(res$table_reaction$Path)))
    expect_true(all(apply(res$table_stat[, 1:2], 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_node, 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_edge[, 1:2], 2, function(x) all(!is.na(x)))))
}

#### Test ####
test_that("Test under normal conditions", {
    for(show_sp in c('all', 'sigClass', 'none')){
        for(show_all_reactions in c(TRUE, FALSE)){
            for(organism in c('human', 'mouse')){
                for(sp_significant in c('pval', 'padj')){
                    for(class_significant in c('pval', 'padj')){
                        res <- nw_lipid_reaction(
                            deSp_se=deSp_se_twoGroup, organism, show_sp,
                            show_all_reactions, sp_significant, sp_p_cutoff=0.05,
                            sp_FC_cutoff=1, class_significant, class_p_cutoff=0.05,
                            class_FC_cutoff=1)
                        lipid_reaction_test(res)
                    }
                }
            }
        }
    }
})

test_that("Test without significant lipids and classes", {
    # res <- nw_lipid_reaction(
    #     deSp_se=deSp_se_twoGroup, organism='mouse', show_sp='sigClass',
    #     show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
    #     sp_FC_cutoff=50, class_significant='pval', class_p_cutoff=0.05,
    #     class_FC_cutoff=50)
    expect_warning(nw_lipid_reaction(
        deSp_se=deSp_se_twoGroup, organism='mouse', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=50, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=50), 'No lipid classes meet these significant cutoffs; please adjust the thresholds you have set.')
    #expect_identical(res, )
})

test_that("Test without significant lipids", {
    res <- nw_lipid_reaction(
        deSp_se=deSp_se_twoGroup, organism='mouse', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=1000, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=1)
    lipid_reaction_test(res)
})

test_that("Test with positive significant lipid species", {
    S4Vectors::metadata(deSp_se_twoGroup)$all_deSp_result %<>%
        dplyr::mutate(log2FC=abs(log2FC))
    res <- nw_lipid_reaction(
        deSp_se=deSp_se_twoGroup, organism='mouse', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=1)
    lipid_reaction_test(res)
})

test_that("Test with negative significant lipid species", {
    S4Vectors::metadata(deSp_se_twoGroup)$all_deSp_result %<>%
        dplyr::mutate(log2FC=-abs(log2FC))
    res <- nw_lipid_reaction(
        deSp_se=deSp_se_twoGroup, organism='mouse', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=1)
    lipid_reaction_test(res)
})


data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_multi <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')

test_that("Test for incorrect input and error message", {
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='humna', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The 'organism' parameter must be either 'human' or 'mouse'."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sig',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The 'show_sp' parameter must be one of 'all', 'sigClass' or 'none'."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='Pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The 'sp_significant' parameter must be either 'pval' or 'padj'."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='Pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The 'class_significant' parameter must be either 'pval' or 'padj'."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=1.1,
            sp_FC_cutoff=2, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The 'sp_p_cutoff' parameter must be a numeric value between 0 and 1."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=1.5,
            class_FC_cutoff=1),
        "The 'class_p_cutoff' parameter must be a numeric value between 0 and 1."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=0, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=1),
        "The minimum value for the 'sp_FC_cutoff' parameter is 1."
    )
    expect_error(
        nw_lipid_reaction(
            deSp_se=deSp_se_twoGroup, organism='human', show_sp='sigClass',
            show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
            sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
            class_FC_cutoff=0),
        "The minimum value for the 'class_FC_cutoff' parameter is 1."
    )
    expect_error(nw_lipid_reaction(
        deSp_se_multi, organism='human', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=1),
        "This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis."
    )
    sub <- processed_se_twoGroup[, SummarizedExperiment::colData(processed_se_twoGroup)$sample_name %in% c("control_01", "hfref_patient_01")]
    deSp_se <- suppressWarnings(
        deSp_twoGroup(
            sub, ref_group='ctrl', test='Wilcoxon test', significant='padj',
            p_cutoff=0.05, FC_cutoff=2, transform='log10'))
    expect_error(nw_lipid_reaction(
        deSp_se, organism='human', show_sp='sigClass',
        show_all_reactions=FALSE, sp_significant='pval', sp_p_cutoff=0.05,
        sp_FC_cutoff=1, class_significant='pval', class_p_cutoff=0.05,
        class_FC_cutoff=1),
    "This function requires each group to have at least two samples."
    )
})
