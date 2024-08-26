library(testthat)
#library(SummarizedExperiment)

data("se_multiGroup")
processed_se <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
post_hoc_p_cutoff <- 0.05

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


