library(testthat)
#library(SummarizedExperiment)
#library(tidyverse)

data("abundance_twoGroup")
data("group_info_twoGroup")
data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_twoGroup <- deSp_twoGroup(
    processed_se_twoGroup, ref_group='ctrl', test='t-test',
    significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

#### Function ####
pathway_activity_test <- function(res){
    # Verify the class of the return objects.
    expect_type(res, 'list')
    expect_s3_class(res$table_edge, 'data.frame')
    expect_s3_class(res$table_node, 'data.frame')
    expect_s3_class(res$table_pathway_score, 'data.frame')
    expect_s3_class(res$table_zScore, 'data.frame')
    # Confirm the existence of required columns.
    expect_true(all(c('from', 'to', 'label', 'width', 'arrows', 'dashes', 'title',
                      'smooth', 'shadow', 'color', 'font.size')
                    %in% colnames(res$table_edge)))
    expect_true(all(c('id', 'label', 'size', 'shape', 'title', 'color', 'border',
                      'font.size', 'shadow') %in% colnames(res$table_node)))
    expect_true(all(c('path_name', 'reaction_chain', 'pathway_score', 'genes')
                    %in% colnames(res$table_pathway_score)))
    expect_true(all(c('from', 'to', 'p', 'z', 'Reactions', 'EC.numbers',
                      'Crossreferences', 'gene') %in% colnames(res$table_zScore)))
    # Confirm whether the data frame consists of numeric values.
    expect_true(all(apply(res$table_zScore[,3:4], 2, function(x) all(is.numeric(x)))))
    expect_true(all(apply(res$table_pathway_score[, 3], 2, function(x) all(is.numeric(x)))))
    # Check NA values
    expect_true(all(apply(res$table_zScore[, 1:2], 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_pathway_score[, 1:2], 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_node, 2, function(x) all(!is.na(x)))))
    expect_true(all(apply(res$table_edge, 2, function(x) all(!is.na(x)))))
}

#### Test ####
test_that("Test under normal conditions", {
    for(organism in c('human', 'mouse')){
        res <- nw_pathway_activity(deSp_se=deSp_se_twoGroup, organism)
        pathway_activity_test(res)
    }
})

test_that("Test under paired samples", {
    sample.list <- c('control_01', 'control_02', 'control_03', 'control_04',
                     'hfref_patient_01', 'hfref_patient_02', 'hfref_patient_03',
                     'hfref_patient_04')
    abundance <- abundance_twoGroup %>%
        dplyr::select(feature, dplyr::all_of(sample.list))
    group_info <- group_info_twoGroup %>%
        dplyr::filter(sample_name %in% sample.list) %>%
        dplyr::arrange(match(sample_name, colnames(abundance))) %>%
        dplyr::mutate(pair=rep(1:4, 2))
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance$feature)
    recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
    abundance %<>% dplyr::filter(feature %in% recognized_lipid)
    goslin_annotation <- parse_lipid %>%
        dplyr::filter(Original.Name %in% recognized_lipid)
    SE_data <- suppressWarnings(
        as_summarized_experiment(
            abundance, goslin_annotation, group_info, n_group='two', paired_sample=TRUE))
    processed_se <- data_process(
        se=SE_data, exclude_missing=TRUE, exclude_missing_pct=70,
        replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage')
    deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
                             significant='pval', p_cutoff=0.05, FC_cutoff=1,
                             transform='log10')
    for(organism in c('human', 'mouse')){
        res <- nw_pathway_activity(deSp_se, organism)
        pathway_activity_test(res)
    }
})


data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
deSp_se_multi <- deSp_multiGroup(
    processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
    significant='pval', p_cutoff=0.05, transform='log10')

test_that("Test for incorrect input data.", {
    expect_error(
        nw_pathway_activity(deSp_se_multi, organism="human"),
        'This function is used for two groups. If your data consists of multiple groups, please subset two of those groups for analysis.'
    )
    expect_error(
        nw_pathway_activity(deSp_se_twoGroup, organism="humna"),
        "The 'organism' parameter must be either 'human' or 'mouse'."
    )
    sub <- processed_se_twoGroup[, SummarizedExperiment::colData(processed_se_twoGroup)$sample_name %in% c("control_01", "hfref_patient_01")]
    deSp_se <- deSp_twoGroup(
        sub, ref_group='ctrl', test='Wilcoxon test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    expect_error(suppressWarnings(
        nw_pathway_activity(deSp_se, organism="human") ),
        "Comparing the means of two independent groups requires each group to have at least two samples."
    )
})
