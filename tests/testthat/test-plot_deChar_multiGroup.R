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
                deChar_se <- deChar_multiGroup(
                    processed_se, char='Total.C', ref_group='ctrl', post_hoc,
                    post_hoc_sig, post_hoc_p_cutoff, transform)
                plot_deChar_se <- suppressWarnings(
                    plot_deChar_multiGroup(deChar_se))
                # Verify the class of objects
                expect_true(is.data.frame(plot_deChar_se$table_barPlot))
                expect_true(is.data.frame(plot_deChar_se$table_linePlot))
                expect_true(is.data.frame(plot_deChar_se$table_boxPlot))
                expect_true(is.data.frame(plot_deChar_se$table_char_index))
                expect_true(is.data.frame(plot_deChar_se$table_index_stat))
                expect_s3_class(plot_deChar_se$static_barPlot, "ggplot")
                expect_s3_class(plot_deChar_se$static_barPlot_sqrt, "ggplot")
                expect_s3_class(plot_deChar_se$static_linePlot, "ggplot")
                expect_s3_class(plot_deChar_se$static_linePlot_sqrt, "ggplot")
                expect_s3_class(plot_deChar_se$static_boxPlot, "ggplot")
                expect_s3_class(plot_deChar_se$interactive_barPlot, "plotly")
                expect_s3_class(plot_deChar_se$interactive_barPlot_sqrt, "plotly")
                expect_s3_class(plot_deChar_se$interactive_linePlot, "plotly")
                expect_s3_class(plot_deChar_se$interactive_linePlot_sqrt, "plotly")
                # Confirm the existence of required columns.
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_deChar_se$table_barPlot)))
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_deChar_se$table_linePlot)))
                expect_true(all(c(
                    'feature', 'sample_name', 'abund', 'label_name', 'original_group_name',
                    'group') %in% colnames(plot_deChar_se$table_boxPlot)))
                expect_true('feature' %in% colnames(plot_deChar_se$table_char_index))
                expect_true(all(c(
                    'feature', 'method', 'statistic', 'pval', 'negLog10pval', 'padj',
                    'negLog10padj') %in% colnames(plot_deChar_se$table_index_stat)))
            }
        }
    }
})

test_that("Test under normal conditions: category characteristics", {
    for(post_hoc in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(post_hoc_sig in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                deChar_se <- deChar_multiGroup(
                    processed_se, char='class', ref_group='ctrl', post_hoc,
                    post_hoc_sig, post_hoc_p_cutoff, transform)
                plot_deChar_se <- plot_deChar_multiGroup(deChar_se)
                expect_message(plot_deChar_multiGroup(deChar_se),
                               info="The 'char' parameter is not a numeric characteristic, so it is unable to perform line plot and box plot.")
                # Verify the class of objects
                expect_true(is.data.frame(plot_deChar_se$table_barPlot))
                expect_true(is.null(plot_deChar_se$table_linePlot))
                expect_true(is.null(plot_deChar_se$table_boxPlot))
                expect_true(is.null(plot_deChar_se$table_char_index))
                expect_true(is.null(plot_deChar_se$table_index_stat))
                expect_s3_class(plot_deChar_se$static_barPlot, "ggplot")
                expect_s3_class(plot_deChar_se$static_barPlot_sqrt, "ggplot")
                expect_true(is.null(plot_deChar_se$static_linePlot))
                expect_true(is.null(plot_deChar_se$static_linePlot_sqrt))
                expect_true(is.null(plot_deChar_se$static_boxPlot))
                expect_s3_class(plot_deChar_se$interactive_barPlot, "plotly")
                expect_s3_class(plot_deChar_se$interactive_barPlot_sqrt, "plotly")
                expect_true(is.null(plot_deChar_se$interactive_linePlot))
                expect_true(is.null(plot_deChar_se$interactive_linePlot_sqrt))
                # Confirm the existence of required columns.
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_deChar_se$table_barPlot)))
            }
        }
    }
})

test_that("Test under normal conditions: ratio characteristics", {
    for(post_hoc in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(post_hoc_sig in c('pval', 'padj')){
            for(transform in c('none', 'log2')){
                deChar_se <- deChar_multiGroup(
                    processed_se, char='Chains Ether/Ester linked ratio',
                    ref_group='ctrl', post_hoc, post_hoc_sig, post_hoc_p_cutoff, transform)
                plot_deChar_se <- plot_deChar_multiGroup(deChar_se)
                expect_message(plot_deChar_multiGroup(deChar_se),
                               info="The 'char' parameter is not a numeric characteristic, so it is unable to perform line plot and box plot.")
                # Verify the class of objects
                expect_true(is.data.frame(plot_deChar_se$table_barPlot))
                expect_true(is.null(plot_deChar_se$table_linePlot))
                expect_true(is.null(plot_deChar_se$table_boxPlot))
                expect_true(is.null(plot_deChar_se$table_char_index))
                expect_true(is.null(plot_deChar_se$table_index_stat))
                expect_s3_class(plot_deChar_se$static_barPlot, "ggplot")
                expect_s3_class(plot_deChar_se$static_barPlot_sqrt, "ggplot")
                expect_true(is.null(plot_deChar_se$static_linePlot))
                expect_true(is.null(plot_deChar_se$static_linePlot_sqrt))
                expect_true(is.null(plot_deChar_se$static_boxPlot))
                expect_s3_class(plot_deChar_se$interactive_barPlot, "plotly")
                expect_s3_class(plot_deChar_se$interactive_barPlot_sqrt, "plotly")
                expect_true(is.null(plot_deChar_se$interactive_linePlot))
                expect_true(is.null(plot_deChar_se$interactive_linePlot_sqrt))
                # Confirm the existence of required columns.
                expect_true(all(c(
                    'characteristic', 'feature', 'significance', 'p', 'group',
                    'mean', 'sd', 'max_error_bar', 'pvalue_text', 'hover',
                    'hovertext') %in% colnames(plot_deChar_se$table_barPlot)))
            }
        }
    }
})

test_that("Test for incorrect input & error message", {
    data("de_data_twoGroup")
    processed_se_two <- data_process(
        de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
        replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage')
    deChar_se_two <- deChar_twoGroup(
        processed_se_two, char="Total.C", ref_group="ctrl", test='t-test',
        significant="padj", p_cutoff=0.05, FC_cutoff=1, transform='log10')
    expect_error(
        plot_deChar_multiGroup(deChar_se_two),
        "This function is used for multiple groups. If your data consists of two groups, please use plot_deChar_twoGroup for analysis."
    )

})
