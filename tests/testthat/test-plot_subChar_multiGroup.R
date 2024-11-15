library(testthat)
#library(SummarizedExperiment)

data("se_multiGroup")
processed_se <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
post_hoc_p_cutoff <- 0.05

data("de_data_twoGroup")
processed_se_two <- data_process(
    de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')

#### Test ####
test_that("Test under normal conditions: continues characteristics", {
    for(post_hoc in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(post_hoc_sig in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                subChar_se <- subChar_multiGroup(
                    processed_se, char='Total.C', subChar='class', ref_group='ctrl',
                    post_hoc, post_hoc_sig, post_hoc_p_cutoff, transform)
                plot_subChar_se <- suppressWarnings(
                    plot_subChar_multiGroup(subChar_se, subChar_feature='CL'))
                # Verify the class of objects
                expect_true(is.data.frame(plot_subChar_se$table_barPlot))
                expect_true(is.data.frame(plot_subChar_se$table_linePlot))
                expect_true(is.data.frame(plot_subChar_se$table_boxPlot))
                expect_true(is.data.frame(plot_subChar_se$table_char_index))
                expect_true(is.data.frame(plot_subChar_se$table_index_stat))
                expect_s3_class(plot_subChar_se$static_barPlot, "ggplot")
                expect_s3_class(plot_subChar_se$static_barPlot_sqrt, "ggplot")
                expect_s3_class(plot_subChar_se$static_linePlot, "ggplot")
                expect_s3_class(plot_subChar_se$static_linePlot_sqrt, "ggplot")
                expect_s3_class(plot_subChar_se$static_boxPlot, "ggplot")
                expect_s3_class(plot_subChar_se$interactive_barPlot, "plotly")
                expect_s3_class(plot_subChar_se$interactive_barPlot_sqrt, "plotly")
                expect_s3_class(plot_subChar_se$interactive_linePlot, "plotly")
                expect_s3_class(plot_subChar_se$interactive_linePlot_sqrt, "plotly")
                expect_true(is.data.frame(plot_subChar_se$table_linePlot))
                expect_true(is.data.frame(plot_subChar_se$table_boxPlot))
                expect_true(is.data.frame(plot_subChar_se$table_char_index))
                expect_true(is.data.frame(plot_subChar_se$table_index_stat))
                # Confirm the existence of required columns.
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_subChar_se$table_barPlot)))
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_subChar_se$table_linePlot)))
                expect_true(all(c(
                    'feature', 'sample_name', 'abund', 'label_name', 'original_group_name',
                    'group') %in% colnames(plot_subChar_se$table_boxPlot)))
                expect_true('feature' %in% colnames(plot_subChar_se$table_char_index))
                expect_true(all(c(
                    'feature', 'method', 'statistic', 'pval', 'negLog10pval', 'padj',
                    'negLog10padj') %in% colnames(plot_subChar_se$table_index_stat)))
            }
        }
    }
})

test_that("plot_subChar_multiGroup handles errors correctly", {
    subChar_se <- subChar_multiGroup(
        processed_se, char='Total.C', subChar='class', ref_group='ctrl',
        post_hoc='One-way ANOVA', post_hoc_sig='pval',
        post_hoc_p_cutoff=post_hoc_p_cutoff, transform='none')

    # Test with invalid subChar_feature
    expect_error(
        plot_subChar_multiGroup(subChar_se, subChar_feature="lipid"),
        "The 'subChar_feature' parameter must be one of the valid subChar features."
    )

    # Test with missing metadata
    invalid_se <- subChar_se
    S4Vectors::metadata(invalid_se)$all_deChar_result <- NULL
    expect_error(
        plot_subChar_multiGroup(invalid_se, subChar_feature="CL"),
        "Correct SummarizedExperiment metadata not found, please use ouput from upstream analysis function."
    )

    # Test with two-group data
    # subChar_se <- subChar_twoGroup(
    #     processed_se_two, char="Total.C", subChar="class",
    #     ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
    #     FC_cutoff=1, transform='log10')
    # expect_error(
    #     plot_subChar_multiGroup(subChar_se, subChar_feature="CL"),
    #     "This function is used for multiple groups. If your data consists of two groups, please use plot_subChar_twoGroup function for analysis."
    # )
})
