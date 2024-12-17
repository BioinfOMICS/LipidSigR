#' @title corr_cor_heatmap
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
#' @param side_color_char Character. A lipid characteristic used for plotting
#' the side color of heatmap. It must be selected from the common list  returned
#' by \code{\link{list_lipid_char}}.
#' @param correlation Character. The method for computing correlation coefficient.
#' Allowed methods include "pearson", "kendall", and "spearman".
#' Default is \code{'pearson'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{1}.
#' @param adjust_p_method Character. The correction method of p-value. Allowed
#' methods include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' and "none". Default is \code{'BH'}.
#' @param cor_coef_cutoff Numeric. Significance of the correlation coefficient.
#' Default is \code{0}.
#' @param distfun Character. The distance measure for computing correlation
#' coefficient (or covariance). Allowed methods include "pearson", "kendall",
#' "spearman". Default is \code{'spearman'}.
#' @param hclustfun Character. The agglomeration method. This should be
#' (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). Default is \code{'average'}.
#' @param heatmap_col Character. The value for clustering. Allow method are
#' "cor_coef" and "statistic". Default is \code{'statistic'}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @param type Character. Specifies the correlation type: 'Sp' for lipid species
#' correlation and 'Char' for lipid characteristic correlation.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("corr_data")
#' processed_se <- data_process(
#'     corr_data, exclude_missing=TRUE, exclude_missing_pct=70, replace_na_method='min',
#'     replace_na_method_ref=0.5, normalization='Percentage')
#' result <- corr_cor_heatmap(processed_se, char=NULL,
#'     condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
#'     side_color_char=NULL, correlation='pearson', significant='pval',
#'     p_cutoff=1, adjust_p_method='BH', cor_coef_cutoff=0, distfun='spearman',
#'     hclustfun='average', heatmap_col='statistic', transform='log10', type='Sp')

corr_cor_heatmap <- function(
        processed_se, char=NULL, condition_col, side_color_char=NULL,
        correlation=c('pearson', 'kendall', 'spearman'),
        significant=c('padj', 'pval'), p_cutoff=1,
        adjust_p_method=c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'),
        cor_coef_cutoff=0, distfun=c('pearson', 'spearman', 'kendall'),
        hclustfun=c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'),
        heatmap_col=c('cor_coef', 'statistic'),
        transform=c('none', 'log10', 'square', 'cube'), type=c('Sp', 'Char')){

    .check_inputSE(processed_se, metadata_list=NULL)
    ## check for condition_col & adjusted_col

    ######
    if (is.null(correlation) | isFALSE(correlation %in% c("pearson", "spearman", "kendall")) ) {
        stop("correlation must be one of 'pearson', 'spearman', or 'kendall'.")
    }
    if (is.null(significant) | isFALSE(significant %in% c('pval', 'padj')) ) {
        stop("significant must be one of 'pval' or 'padj'.")
    }
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    if (is.null(adjust_p_method) | isFALSE(adjust_p_method %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')) ) {
        stop("adjust_p_method must be one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'.")
    }
    if (!is.numeric(cor_coef_cutoff) | isFALSE(.check_numeric_range(cor_coef_cutoff, 0, 1)) ) {
        stop("cor_coef_cutoff must be a numeric value between 0 and 1.")
    }
    if (is.null(distfun) | isFALSE(distfun %in% c('pearson', 'spearman', 'kendall')) ) {
        stop("distfun must be one of 'pearson', 'spearman', or 'kendall'.")
    }
    if (is.null(hclustfun) | isFALSE(hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')) ) {
        stop("hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    }
    if (is.null(heatmap_col) | isFALSE(heatmap_col %in% c('cor_coef', 'statistic')) ) {
        stop("heatmap_col must be 'cor_coef' or 'statistic'.")
    }
    if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
        stop("transform must be one of 'none', 'log10', 'square', or 'cube'.")
    }
    ## check for type
    if (!is.null(type) && type=="Sp") {
        if (!is.null(char)) {
            stop("If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.")
        }
        if (isFALSE(.check_char(processed_se, side_color_char, type='common'))) {
            stop("Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.")
        }
        abundance_raw <- .extract_df(processed_se, type = "abundance")
        lipid_char_table <- .extract_df(processed_se,type='lipid')
    } else if (!is.null(type) && type=="Char") {
        if (!is.null(side_color_char)) {
            stop("If type is set to 'Char' for running lipid characteristics correlation analysis, side_color_char must be set to NULL.")
        }
        char_list <- list_lipid_char(processed_se)$common_list
        if(is.null(char) | isFALSE(char %in% char_list)){
            stop("Wrong char input, you can view the available char list by list_lipid_char function.")
        }
        char_se <- convert_sp2char(processed_se, transform='none')
        char_data <- .extract_df(char_se, type = "abundance")
        lipid_char_table <- .extract_df(char_se,type='lipid')
        abundance_raw <- merge(
            lipid_char_table, char_data, by.x='row.names', by.y='feature') %>%
            dplyr::filter(characteristic == char) %>%
            dplyr::select(-Row.names, -characteristic)
    } else {
        stop("type must be one of 'Sp' or 'Char'.")
    }
    abundance <- .transform(abundance_raw, transform)
    condition_table <- .extract_df(processed_se,type='group') %>%
        dplyr::select(c(sample_name, all_of(condition_col)) ) %>% as.data.frame()

    .check_imputation(abundance)

    CHAR <- colnames(abundance)[1]
    colnames(abundance)[1] <- 'feature'
    abundance %<>% tidyr::gather(-feature, key='sample_name', value='value') %>%
        tidyr::spread(key='feature', value='value') %>% dplyr::arrange(sample_name)

    condition_table <- condition_table %>% dplyr::arrange(sample_name)
    #---------Cor-----------------------------
    cor_table_all <- .corr(
        abundance, condition_table, correlation, p_cutoff, cor_coef_cutoff, adjust_p_method)
    cor_table_sig <- cor_table_all %>%
        dplyr::filter(!!rlang::sym(paste0("sig_", significant))=='yes')
    #---------heatmap-----------------------------
    heatmap <- suppressWarnings(.corr_heatmap(
        cor_table_all, cor_table_sig, heatmap_col, lipid_char_table, side_color_char,
        significant, distfun, hclustfun, cor_type="cor"))
    colnames(cor_table_all)[2] <- CHAR
    colnames(cor_table_sig)[2] <- CHAR

    return(list(
        all_correlation_result=cor_table_all, sig_correlation_result=cor_table_sig,
        interactive_heatmap=heatmap$in.heatmap, static_heatmap=heatmap$static_hm,
        heatmap_matrix=heatmap$plot.mat))
}

.corr <- function(abundance, condition_table, correlation, p_cutoff, cor_coef_cutoff, adjust_p_method) {
    # Perform correlation analysis for each clinical factor
    cor_table_all <- purrr::map_dfr(2:ncol(condition_table), function(a) {
        # Extract clinical factor column
        clin_factor_name <- colnames(condition_table)[a]
        clin_factor_col <- condition_table[[a]]
        # Perform correlation for each feature in abundance
        cor_results <- purrr::map_dfr(2:ncol(abundance), function(b) {
            feature_name <- colnames(abundance)[b]
            feature_col <- abundance[[b]]
            # Perform correlation test with error handling
            cor_test <- tryCatch(
                cor.test(feature_col, clin_factor_col, method=correlation),
                error=function(e){NULL})
            # Extract or set correlation results
            if (!is.null(cor_test)) {
                tibble::tibble(
                    clin_factor=clin_factor_name, feature=feature_name,
                    method=correlation, cor_coef=cor_test$estimate,
                    statistic=cor_test$statistic, pval=cor_test$p.value)
            } else {
                tibble::tibble(
                    clin_factor=clin_factor_name, feature=feature_name,
                    method=correlation, cor_coef=NA, statistic=NA, pval=NA)
            }
        })
        # Add adjusted p-values to the results
        cor_results %>%
            dplyr::mutate(
                padj=p.adjust(pval, method=adjust_p_method, n=length(pval)),
                sig_pval=ifelse(abs(cor_coef)>cor_coef_cutoff & pval<p_cutoff, 'yes','no'),
                sig_padj=ifelse(abs(cor_coef)>cor_coef_cutoff & padj<p_cutoff, 'yes','no')) %>%
            dplyr::mutate(pval=round(pval, 5), padj=round(padj, 5)) %>% as.data.frame()
    })
    return(cor_table_all=cor_table_all)
}
