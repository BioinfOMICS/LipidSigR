library(testthat)
library(SummarizedExperiment)

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage')
p_cutoff <- 0.05

DE_multi_sig_lipid <- function(DE.res, p_cutoff){
    # Verify if the significant lipids indeed represent significance.
    expect_equal(DE.res$feature[which(DE.res$pval < p_cutoff)],
                 DE.res$feature[which(DE.res$sig_pval == 'yes')])
    expect_equal(DE.res$feature[which(DE.res$padj < p_cutoff)],
                 DE.res$feature[which(DE.res$sig_padj == 'yes')])
}

#### Test ####
test_that("Test under normal conditions", {
    for(test in c('One-way ANOVA', 'Kruskal–Wallis test')){
        for(significant in c('pval', 'padj')){
            for(transform in c('none', 'log10', 'square', 'cube')){
                deSp_se <- deSp_multiGroup(
                    processed_se=processed_se_multiGroup, ref_group='ctrl', test,
                    significant, p_cutoff, transform)
                plot_deSp_se <- plot_deSp_multiGroup(deSp_se)
                # Verify the class of objects
                expect_true(is.data.frame(plot_deSp_se$table_de_lipid))
                expect_true(is.data.frame(plot_deSp_se$table_dotPlot))
                expect_s3_class(plot_deSp_se$static_de_lipid, "ggplot")
                expect_s3_class(plot_deSp_se$static_dotPlot, "ggplot")
                expect_s3_class(plot_deSp_se$interactive_de_lipid, "plotly")
                expect_s3_class(plot_deSp_se$interactive_dotPlot, "plotly")
                # Confirm the existence of required columns.
                expect_true(all(c('feature', 'pval', 'negLog10pval', 'padj',
                                  'negLog10padj', 'sig_pval', 'sig_padj', 'hover')
                                %in% colnames(plot_deSp_se$table_de_lipid)))
                expect_false('log2FC' %in% colnames(plot_deSp_se$table_de_lipid))
                expect_true(all(c('feature', 'pval', 'negLog10pval', 'padj', 'class',
                                  'negLog10padj', 'sig_pval', 'sig_padj', 'hover')
                                %in% colnames(plot_deSp_se$table_dotPlot)))
                expect_false('log2FC' %in% colnames(plot_deSp_se$table_dotPlot))
                # Confirm the correctness of significance.
                DE_multi_sig_lipid(plot_deSp_se$table_de_lipid, p_cutoff)
                DE_multi_sig_lipid(plot_deSp_se$table_dotPlot, p_cutoff)
            }
        }
    }
})






