library(testthat)
library(SummarizedExperiment)

data("de_data_twoGroup")
data("se_multiGroup")
p_cutoff <- 0.05
FC_cutoff <- 1

## Two Group
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
## Multi Group
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

#### Function ####
two_char_heatmap_se_format <- function(se) {
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
    # expect_true("metadata" %in% slotNames(se),
    #             info = "metadata slot is missing")
    # Check the structure of assays, rowData, and colData
    # Check expression data exists and is a matrix
    expect_true(is.matrix(assay(se)),
                info = "Expression data is not a matrix")
    expect_true(all(dim(assay(se)) > 0),
                info = "Expression data has non-positive dimensions")
}
heatmap_chain_db_test <- function(result, p_cutoff, n_group){
    # Verify the SE object returned.
    two_char_heatmap_se_format(se=result$chain_db_se)
    # Explicitly force the evaluation of the plot.
    expect_s3_class(result$static_heatmap, "ggplot")
    # Verify if it is a data frame.
    expect_true(is.data.frame(result$table_heatmap))
    expect_true(is.data.frame(result$processed_abundance))
    expect_true(is.data.frame(result$transformed_abundance))
    # Confirm the existence of required columns.
    if(n_group == 'two'){
        expect_true(all(c('feature', 'method', 'FC', 'log2FC', 'pval', 'negLog10pval',
                          'padj', 'negLog10padj', 'sig_pval', 'sig_padj', 'Chain', 'DB')
                        %in% colnames(result$table_heatmap)))
    }else if(n_group == 'multiple'){
        expect_true(all(c('feature', 'method', 'pval', 'negLog10pval', 'padj',
                          'negLog10padj', 'sig_pval', 'sig_padj', 'Chain', 'DB')
                        %in% colnames(result$table_heatmap)))
        expect_false(all(c('FC', 'log2FC') %in% colnames(result$table_heatmap)))
    }
    expect_true('feature' %in% colnames(result$processed_abundance))
    expect_true('feature' %in% colnames(result$transformed_abundance))
    # Confirm the correctness of significance.
    expect_equal(result$table_heatmap$feature[which(result$table_heatmap$sig_pval == 'yes')],
                 result$table_heatmap$feature[which(result$table_heatmap$pval < p_cutoff)])
    expect_equal(result$table_heatmap$feature[which(result$table_heatmap$sig_padj == 'yes')],
                 result$table_heatmap$feature[which(result$table_heatmap$padj < p_cutoff)])
    # Confirm if the lipid features are the same.
    expect_equal(result$table_heatmap$feature, result$processed_abundance$feature)
    expect_equal(result$table_heatmap$feature, result$transformed_abundance$feature)
}

#### Test ####
test_that("Test under normal conditions without 'charFeature': two groups", {
    for(test in c('t-test', 'Wilcoxon test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                res <- heatmap_chain_db(
                    processed_se=processed_se_twoGroup, char='class', char_feature=NULL,
                    ref_group='ctrl', test, significant, p_cutoff, FC_cutoff, transform=transform)
                heatmap_chain_db_test(result=res$total_chain, p_cutoff, n_group='two')
                heatmap_chain_db_test(result=res$each_chain, p_cutoff, n_group='two')
            }
        }
    }
})

test_that("Test under normal conditions without 'charFeature': multi groups", {
    for(test in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                res <- heatmap_chain_db(
                    processed_se=processed_se_multiGroup, char='class', char_feature=NULL,
                    ref_group='ctrl', test, significant, p_cutoff, FC_cutoff, transform=transform)
                heatmap_chain_db_test(result=res$total_chain, p_cutoff, n_group='multiple')
                heatmap_chain_db_test(result=res$each_chain, p_cutoff, n_group='multiple')
            }
        }
    }
})

test_that("Test under normal conditions with 'charFeature': two groups", {
    for(test in c('t-test', 'Wilcoxon test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                res <- heatmap_chain_db(
                    processed_se=processed_se_twoGroup, char='class', char_feature='PC',
                    ref_group='ctrl', test, significant, p_cutoff, FC_cutoff, transform=transform)
                heatmap_chain_db_test(result=res$total_chain, p_cutoff, n_group='two')
                heatmap_chain_db_test(result=res$each_chain, p_cutoff, n_group='two')
            }
        }
    }
})

test_that("Test under normal conditions with 'charFeature': multi groups", {
    for(test in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                res <- heatmap_chain_db(
                    processed_se=processed_se_multiGroup, char='class', char_feature='PC',
                    ref_group='ctrl', test, significant, p_cutoff, FC_cutoff, transform=transform)
                heatmap_chain_db_test(result=res$total_chain, p_cutoff, n_group='multiple')
                heatmap_chain_db_test(result=res$each_chain, p_cutoff, n_group='multiple')
            }
        }
    }
})

test_that("Test special conditions: two groups", {
    processed_se <- processed_se_twoGroup[1:5,]
    res1 <- heatmap_chain_db(
        processed_se, char='class', char_feature=NULL,
        ref_group='ctrl', test='t-test', significant='pval',
        p_cutoff, FC_cutoff, transform='log10')
    heatmap_chain_db_test(result=res1$total_chain, p_cutoff, n_group='two')
    expect_true(is.character(res1$each_chain))
})

test_that("Test error conditions", {
    # Invalid char
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='invalid_char', ref_group='ctrl'))

    # Invalid char_feature
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', char_feature='invalid_feature', ref_group='ctrl'))

    # Invalid ref_group
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='invalid_group'))

    # Invalid test for two groups
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', test='One-way ANOVA'))

    # Invalid test for multi groups
    expect_error(heatmap_chain_db(processed_se_multiGroup, char='class', ref_group='ctrl', test='t-test'))

    # Invalid significant
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', significant='invalid'))

    # Invalid p_cutoff
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', p_cutoff=2))

    # Invalid FC_cutoff
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', FC_cutoff=0))

    # Invalid transform
    expect_error(heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', transform='invalid'))
})

test_that("Test edge cases", {
    # Empty SummarizedExperiment
    empty_se <- SummarizedExperiment(assays=matrix(nrow=0, ncol=0))
    expect_error(
        heatmap_chain_db(empty_se, char='class', ref_group='ctrl'),
        "The input must be a SummarizedExperiment object construct by as_summarized_experiment function or output from upstream analysis function.")

    # Single feature
    single_feature_se <- processed_se_twoGroup[1,]
    expect_error(
        heatmap_chain_db(single_feature_se, char='class', ref_group='ctrl'),
        "Incorrect SummarizedExperiment structure, please construct by as_summarized_experiment function."
    )

    # All p-values above cutoff
    expect_error(
        heatmap_chain_db(processed_se_twoGroup, char='class', ref_group='ctrl', p_cutoff=1e-10),
        "'One-way ANOVA' and 'Kruskal–Wallis test' are for multiple group comparisons. Please choose either 't-test' or 'Wilcoxon test' for two group analysis."
    )
    res_high_p <- suppressWarnings(
        heatmap_chain_db(
            processed_se_twoGroup, char='class', char_feature=NULL, ref_group='ctrl',
            test='t-test', significant='pval', p_cutoff=1e-10,  FC_cutoff=1, transform='log10')
    )
    if (!is.character(res_high_p$total_chain) ) {
        heatmap_chain_db_test(res_high_p$total_chain, p_cutoff=1e-10, n_group="two")
    } else {
        expect_identical(
            res_high_p$total_chain,
            "This lipid characteristic or characteristic feature does not include any lipids for total chain analysis. Therefore, the 'total_chain' object is NULL."
        )
    }
    if (!is.character(res_high_p$each_chain) ) {
        heatmap_chain_db_test(res_high_p$total_chain, p_cutoff=1e-10, n_group="two")
    } else {
        expect_identical(
            res_high_p$each_chain,
            "This lipid characteristic or characteristic feature does not include any lipids for total chain analysis. Therefore, the 'total_chain' object is NULL."
        )
    }


    # # All fold changes below cutoff (for two groups)
    # res_high_fc <- suppressWarnings(
    #     heatmap_chain_db(
    #         processed_se_twoGroup, char='class', char_feature=NULL, ref_group='ctrl',
    #         test='t-test', significant='pval', p_cutoff=0.05,  FC_cutoff=1000, transform='log10')
    # )
    # if (!is.character(res_high_fc$total_chain) ) {
    #     heatmap_chain_db_test(res_high_fc$total_chain, p_cutoff=0.05, n_group="two")
    # } else {
    #     expect_identical(
    #         res_high_fc$total_chain,
    #         "This lipid characteristic or characteristic feature does not include any lipids for total chain analysis. Therefore, the 'total_chain' object is NULL."
    #     )
    # }
    # if (!is.character(res_high_fc$each_chain) ) {
    #     heatmap_chain_db_test(res_high_fc$each_chain, p_cutoff=0.05, n_group="two")
    # } else {
    #     expect_identical(
    #         res_high_fc$each_chain,
    #         "This lipid characteristic or characteristic feature does not include any lipids for total chain analysis. Therefore, the 'total_chain' object is NULL."
    #     )
    # }
})

