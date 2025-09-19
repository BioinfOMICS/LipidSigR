library(testthat)
# library(SummarizedExperiment)
# library(tidyverse)
data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')


#### Function ####
DE_multi_box_test <- function(res, lipid){
    # Verify if it is a list.
    expect_type(res, 'list')
    # Explicitly force the evaluation of the plot.
    expect_s3_class(res$static_boxPlot, "ggplot")
    # Verify if it is a data frame.
    expect_s3_class(res$table_boxplot, "data.frame")
    expect_s3_class(res$table_stat, "data.frame")
    # Confirm the existence of required columns.
    expect_true(all(c('feature', 'sample_name', 'abund', 'label_name',
                      'original_group_name', 'group')
                    %in% colnames(res$table_boxplot)))
    expect_true(all(c('feature', 'contrast', 'method', 'statistic', 'pval',
                      'post_hoc_method')
                    %in% colnames(res$table_stat)))
    # Confirm if the lipid features are the same.
    expect_equal(lipid, unique(res$table_boxplot$feature))
    expect_equal(lipid, unique(res$table_stat$feature))
    # Numeric value.
    expect_true(is.numeric(res$table_boxplot$abund))
    expect_true(all(sapply(res$table_stat[, c(4:5, 7)], is.numeric)))
}


#### Test ####
test_that("Test results from data_process.R", {
    for(test in c('One-way ANOVA', 'Kruskal-Wallis test')){
        for(transform in c('none', 'log10', 'square', 'cube')){
            if(test == 'One-way ANOVA'){
                post_hoc_sig_method <- 'padj'
            }else{
                post_hoc_sig_method <- c('pval', 'padj')
            }
            for(post_hoc_sig in post_hoc_sig_method){
                res <- boxPlot_feature_multiGroup(
                    processed_se=processed_se_multiGroup,
                    feature='PE O- 17:0;0_20:3;0', ref_group='ctrl',
                    test, post_hoc_sig, transform)
                # if(is.list(res)){
                    DE_multi_box_test(res, lipid='PE O- 17:0;0_20:3;0')
                # }else{
                #     expect_identical(res, 'Unable to perform statistical tests for this feature; therefore, results could not be generated.')
                # }
            }
        }
    }
})

test_that("Test results from heatmap_chain_db.R", {
    for(test in c('One-way ANOVA', 'Kruskal-Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                if(test == 'One-way ANOVA'){
                    post_hoc_sig_method <- 'padj'
                }else{
                    post_hoc_sig_method <- c('pval', 'padj')
                }
                for(post_hoc_sig in post_hoc_sig_method){
                    twoChar <- heatmap_chain_db(
                        processed_se=processed_se_multiGroup, char='class',
                        char_feature='PC', ref_group='ctrl', test, significant,
                        p_cutoff=0.05, FC_cutoff=1, transform)
                    ## Total chain
                    res <- boxPlot_feature_multiGroup(
                        processed_se=twoChar$total_chain$chain_db_se,
                        feature='36:4', ref_group='ctrl',
                        test, post_hoc_sig, transform)
                    DE_multi_box_test(res, lipid='36:4')
                    res2 <- boxPlot_feature_multiGroup(
                        processed_se=twoChar$each_chain$chain_db_se,
                        feature='16:0', ref_group='ctrl',
                        test, post_hoc_sig, transform)
                    DE_multi_box_test(res2, lipid='16:0')
                }

            }
        }
    }
})


data("de_data_twoGroup")
processed_se_two <- data_process(
    de_data_twoGroup, exclude_missing=TRUE,
    exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')


test_that("Test that message return for invalid input", {
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_multiGroup, feature='ABC',
            ref_group='ctrl', test="One-way ANOVA", post_hoc_sig='padj', transform='log10'),
        "'feature' must be one of the lipids in your lipid abundance data; please revise your input.")
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_multiGroup, feature='PE O- 17:0;0_20:3;0',
            ref_group='AA1', test="One-way ANOVA", post_hoc_sig='padj', transform='log10'),
        "ref_group must be one of the group names in the 'group' column of the group information table.")
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_multiGroup, feature='PE O- 17:0;0_20:3;0',
            ref_group='ctrl', test="One-way-ANOVA", post_hoc_sig='padj', transform='log10'),
        "The 'test' parameter must be either 'One-way ANOVA' or 'Kruskal-Wallis test'.")
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_multiGroup, feature='PE O- 17:0;0_20:3;0',
            ref_group='ctrl', test="Kruskal-Wallis test", post_hoc_sig='p-val', transform='log10'),
        "The 'post_hoc_sig' parameter must be either 'pval' or 'padj'.")
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_multiGroup, feature='PE O- 17:0;0_20:3;0',
            ref_group='ctrl', test="One-way ANOVA", post_hoc_sig='padj', transform='log2'),
        "The transform parameter must be one of 'none', 'log10', 'square' or 'cube'.")
    expect_error(
        boxPlot_feature_multiGroup(
            processed_se=processed_se_two, feature='PC O- 16:1;0-20:3;0',
            ref_group='ctrl', test="One-way ANOVA", post_hoc_sig='padj', transform='log10'),
        "This function is used for multiple groups. If your data consists of two groups, please use boxPlot_feature_twoGroup\\(\\) for analysis.")
})
