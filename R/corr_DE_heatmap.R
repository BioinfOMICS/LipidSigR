#' @title corr_DE_heatmap
#' @description This function provides a summary by correlation coefficient of
#' the relationship between clinical features and lipid species or characteristics,
#' indicating its strength and whether it is positive or negative.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @param condition_col Character.The column names used to extract the condition
#' table from the group information table, including clinical conditions such as
#' disease status or gene dependency scores.
#' @param adjust_p_method Character. The correction method of p-value. Allowed
#' methods include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' and "none". Default is \code{'BH'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{1}.
#' @param ParametricTest Logical. Whether to use parametric tests for differential
#' expression analysis. If \code{TRUE} (default), uses t-test for two-group comparison
#' and One-way ANOVA for multi-group comparison. If \code{FALSE}, uses Wilcoxon test
#' for two-group comparison and Kruskal-Wallis test for multi-group comparison.
#' @param side_color_char Character. A lipid characteristic used for plotting
#' the side color of heatmap. It must be selected from the common list  returned
#' by \code{\link{list_lipid_char}}.
#' @param distfun Character. The distance measure for computing correlation
#' coefficient (or covariance). Allowed methods include "pearson", "kendall",
#' "spearman". Default is \code{'spearman'}.
#' @param hclustfun Character. The agglomeration method. This should be
#' (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). Default is \code{'average'}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @param type Character. Specifies the correlation type: 'Sp' for lipid species
#' correlation and 'Char' for lipid characteristic correlation. Default is \code{'Sp'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("corr_data")
#' processed_se <- data_process(
#'     corr_data, exclude_missing=TRUE, exclude_missing_pct=70, replace_na_method='min',
#'     replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
#' result <- corr_DE_heatmap(
#'     processed_se, char=NULL, condition_col=c('Sex', 'Smoking', 'Exacerbations'),
#'     adjust_p_method='BH', significant='pval', p_cutoff=1, side_color_char=NULL,
#'     distfun='spearman', hclustfun='average', transform='log10', type='Sp')


corr_DE_heatmap <- function(
        processed_se, char=NULL, condition_col, side_color_char=NULL,
        significant='pval', p_cutoff=1, ParametricTest=TRUE,
        adjust_p_method='BH', distfun='spearman', hclustfun='average',
        transform='log10', type='Sp'){

    .check_inputSE(processed_se, metadata_list=NULL)
    ## check for condition_col & adjusted_col
    group_info <- .extract_df(processed_se, type='group')
    .check_df_corrCol(group_info, condition_col, adjusted_col=NULL, cor_type='lr')
    match.arg(significant, c('pval', 'padj'))
    match.arg(adjust_p_method, c('BH', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BY', 'fdr', 'none'))
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    match.arg(distfun, c('spearman', 'pearson', 'kendall'))
    match.arg(hclustfun, c('average', 'ward.D', 'ward.D2', 'single', 'complete', 'mcquitty', 'median', 'centroid'))
    match.arg(transform, c('log10', 'square', 'cube', 'none'))
    match.arg(type, c('Sp', 'Char'))

    if (type == "Sp") {
        if (!is.null(char)) {
            stop("If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.")
        }
        if (isFALSE(.check_char(processed_se, side_color_char, type='common'))) {
            stop("Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.")
        }
        abundance <- .extract_df(processed_se, type = "abundance")
        lipid_char_table <- .extract_df(processed_se,type='lipid')
    } else if (type == "Char") {
        if (!is.null(side_color_char)) {
            stop("If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL.")
        }
        char_list <- list_lipid_char(processed_se)$common_list
        if(is.null(char) | isFALSE(char %in% char_list)){
            stop("Wrong char input, you can view the available char list by list_lipid_char function.")
        }
        char_se <- convert_sp2char(processed_se, transform='none')
        char_data <- .extract_df(char_se, type = "abundance")
        lipid_char_table <- .extract_df(char_se, type='lipid')
        abundance <- merge(
            lipid_char_table, char_data, by.x='row.names', by.y='feature') %>%
            dplyr::filter(characteristic == char) %>%
            dplyr::select(-Row.names, -characteristic)
    }
    conCol <- c("sample_name", condition_col)
    condition_table <- .extract_df(processed_se,type='group') %>%
        dplyr::select(all_of(conCol)) %>% as.data.frame()
    .check_imputation(abundance)

    uniqueCounts <- sapply(condition_table[-1], function(col) length(unique(col)))
    DE_table_all <- list(length = (ncol(condition_table)-1))
    DE_table_sig <- list(length = (ncol(condition_table)-1))
    clin_var <- vector(length = (ncol(condition_table)-1))
    for(a in 2:ncol(condition_table)){

        k <- condition_table[c(1,a)]
        k <- k %>% dplyr::mutate(label_name=sample_name, pair=NA)
        clin_var[a-1] <- colnames(k)[2]
        colnames(k)[2] <- 'group'
        if(length(unique(stats::na.omit(k$group))) == 2){
            group_levels <- sort(unique(stats::na.omit(k$group)))
            k$group <- ifelse(k$group == group_levels[1], 'ctrl',
                              ifelse(k$group == group_levels[2], 'exp', NA))
        }

        non_na_group <- k[(!is.na(k$group)),] %>% .$sample_name
        k1 <- k %>% dplyr::filter(sample_name %in% non_na_group)
        non_na_group2 <- colnames(abundance) %in% non_na_group
        non_na_group2[1] <- TRUE
        abundance2 <- abundance[non_na_group2]
        paired_sample <- FALSE


        if(length(unique(k$group)) == 2){
            abundance_tab <- .mean_sd_twoGroup(abundance2, group_info=k1)
            ## Data transformation
            abundance_cli <- .transform(abundance2, transform)
            stat_table <- .stat_twoGroup(
                abundance_cli, group_info=k1, paired_sample,
                test=ifelse(ParametricTest, 't-test', 'Wilcoxon test'))
            ## result table
            res_table <- abundance_tab %>%
                dplyr::left_join(stat_table, by='feature') %>%
                as.data.frame()
            diff_table <- .sig_feature(res_table, significant, p_cutoff, FC_cutoff=NULL)
            DE_table_all[[a-1]] <- diff_table$all_table %>%
                dplyr::select(-FC, -log2FC) %>% dplyr::mutate(clin_factor=clin_var[a-1]) %>%
                dplyr::select(-dplyr::starts_with("mean_"), -dplyr::starts_with("sd_"))
            DE_table_sig[[a-1]] <- diff_table$sig_table %>%
                dplyr::select(-FC, -log2FC) %>% dplyr::mutate(clin_factor=clin_var[a-1])
        }else{
            abund.mean <- .mean_sd_multiGroup(abundance2, group_info=k1)
            ## Data transformation
            abundance_cli <- .transform(abundance2, transform)
            comp.res <- .stat_multiGroup(
                abundance_cli, group_info=k1,
                test=ifelse(ParametricTest, 'One-way ANOVA', 'Kruskal-Wallis test'))
            stat_table <- abund.mean %>% dplyr::left_join(comp.res, by='feature') %>%
                as.data.frame()
            res.tab <- .sig_feature(stat_table, significant, p_cutoff)
            DE_table_all[[a-1]] <- res.tab$all_table %>% dplyr::mutate(clin_factor=clin_var[a-1]) %>%
                dplyr::select(-dplyr::starts_with("mean_"), -dplyr::starts_with("sd_"))
            DE_table_sig[[a-1]] <- res.tab$sig_table %>% dplyr::mutate(clin_factor=clin_var[a-1])
        }
    }

    keep_cols <- c("feature", "method", "statistic", "pval", "negLog10pval", "padj", "negLog10padj", "sig_pval", "sig_padj", "clin_factor")

    DE_table_all <- DE_table_all %>%
        lapply(function(df) df[, intersect(keep_cols, colnames(df)), drop = FALSE]) %>%
        dplyr::bind_rows()

    DE_table_sig <- DE_table_sig %>%
        lapply(function(df) df[, intersect(keep_cols, colnames(df)), drop = FALSE]) %>%
        dplyr::bind_rows()

    heatmap <- suppressWarnings(.corr_heatmap(
        DE_table_all, DE_table_sig, heatmap_col='statistic', lipid_char_table, side_color_char,
        significant, distfun, hclustfun, cor_type="lr"))

    if(!is.null(char)){
        DE_table_all %<>% dplyr::rename(!!char := feature)
        DE_table_sig %<>% dplyr::rename(!!char := feature)
    }

    return(
        list(
            all_correlation_result=DE_table_all, sig_correlation_result=DE_table_sig,
            interactive_heatmap=heatmap$in.heatmap, static_heatmap=heatmap$static_hm,
            heatmap_matrix=heatmap$plot.mat))
}
