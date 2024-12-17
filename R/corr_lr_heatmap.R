#' @title corr_lr_heatmap
#' @description This function uses multiple explanatory variables to predict the
#' outcome of a continuous response variable using linear regression. It enables
#' researchers to estimate associations between lipid levels and clinical features.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @param condition_col Character.The column names used to extract the condition
#' table from the group information table, including clinical conditions such as
#' disease status or gene dependency scores.
#' @param adjusted_col Character. The column names used to extract the adjusted table from
#' the group information table, including additional variables to be incorporated
#' into the algorithm for adjusting confounding effects.
#' @param side_color_char Character. A lipid characteristic used for plotting
#' the side color of heatmap. It must be selected from the common list  returned
#' by \code{\link{list_lipid_char}}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{1}.
#' @param adjust_p_method Character. The correction method of p-value. Allowed
#' methods include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' and "none". Default is \code{'BH'}.
#' @param distfun Character. The distance measure for computing correlation
#' coefficient (or covariance). Allowed methods include "pearson", "kendall",
#' "spearman". Default is \code{'spearman'}.
#' @param hclustfun Character. The agglomeration method. This should be
#' (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). Default is \code{'centroid'}.
#' @param heatmap_col Character. The value for clustering. Allow method are
#' "beta_coef" and "t_statistic". Default is \code{'t_statistic'}.
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
#' result <- corr_lr_heatmap(processed_se, char=NULL,
#'     condition_col=c("FEV1_FVC", "Emphysema", "Exacerbations"),
#'     adjusted_col=c("Age", "Sex", "Smoking", "BMI", "FEV1"),
#'     side_color_char=NULL, significant='pval', p_cutoff=0.05,
#'     adjust_p_method='BH', distfun='spearman', hclustfun='centroid',
#'     heatmap_col='t_statistic', transform='log10', type='Sp')

corr_lr_heatmap <- function(
        processed_se, char=NULL, condition_col, adjusted_col, side_color_char,
        significant=c('padj', 'pval'), p_cutoff=0.05,
        adjust_p_method=c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'),
        distfun=c('pearson', 'spearman', 'kendall'),
        hclustfun=c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'),
        heatmap_col=c('beta_coef', 't_statistic'),
        transform=c('none', 'log10', 'square', 'cube'), type=c('Sp', 'Char')){

    .check_inputSE(processed_se, metadata_list=NULL)
    ## check for condition_col & adjusted_col

    ######
    if (is.null(significant) | isFALSE(significant %in% c('pval', 'padj')) ) {
        stop("significant must be one of 'pval' or 'padj'.")
    }
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    if (is.null(adjust_p_method) | isFALSE(adjust_p_method %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')) ) {
        stop("adjust_p_method must be one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'.")
    }
    if (is.null(distfun) | isFALSE(distfun %in% c('pearson', 'spearman', 'kendall')) ) {
        stop("distfun must be one of 'pearson', 'spearman', or 'kendall'.")
    }
    if (is.null(hclustfun) | isFALSE(hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')) ) {
        stop("hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    }
    if (is.null(heatmap_col) | isFALSE(heatmap_col %in% c('beta_coef', 't_statistic')) ) {
        stop("heatmap_col must be 'beta_coef' or 't_statistic'.")
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
    if (!is.null(adjusted_col)) {
        adjusted_table <- .extract_df(processed_se,type='group') %>%
            dplyr::select(c(sample_name, all_of(adjusted_col)) ) %>% as.data.frame()
    } else {
        adjusted_table <- NULL
    }
    .check_imputation(abundance)

    CHAR <- colnames(abundance)[1]
    colnames(abundance)[1] <- 'feature'
    abundance <- abundance %>% tidyr::gather(-feature, key='sample_name', value='value') %>%
        tidyr::spread(key='feature', value='value') %>% dplyr::arrange(sample_name)

    condition_table <- condition_table %>% dplyr::arrange(sample_name)
    if(!is.null(adjusted_table)){
        adjusted_table <- adjusted_table %>% dplyr::arrange(sample_name)
    }
    abundance[-1] <- scale(abundance[-1]) %>% as.data.frame()

    lr_table_all <- .lr(abundance, condition_table, adjusted_table, p_cutoff, adjust_p_method)
    lr_table_sig <- lr_table_all %>%
        dplyr::filter(!!rlang::sym(paste0("sig_", significant))=='yes')
    #---------heatmap-----------------------------
    heatmap <- suppressWarnings(.corr_heatmap(
        lr_table_all, lr_table_sig, heatmap_col, lipid_char_table, side_color_char,
        significant, distfun, hclustfun, cor_type="lr"))

    colnames(lr_table_all)[2] <- CHAR
    colnames(lr_table_sig)[2] <- CHAR

    return(list(
        all_correlation_result=lr_table_all, sig_correlation_result=lr_table_sig,
        interactive_heatmap=heatmap$in.heatmap, static_heatmap=heatmap$static_hm,
        heatmap_matrix=heatmap$plot.mat))
}


.lr <- function(abundance, condition_table, adjusted_table, p_cutoff, adjust_p_method) {
    regression_for_factor <- function(factor_index) {
        # Prepare clinical factor column
        clin_factor_name <- colnames(condition_table)[factor_index]
        clin_factor_col <- condition_table[factor_index]
        # Perform regression for each feature
        regression_results <- purrr::map_dfr(2:ncol(abundance), function(a) {
            # Prepare data for regression
            if (!is.null(adjusted_table)) {
                data_lm <- cbind(clin_factor_col, abundance[a], adjusted_table[-1])
            } else {
                data_lm <- cbind(clin_factor_col, abundance[a])
            }
            colnames(data_lm)[1] <- 'clin_var'

            # Perform regression with error handling
            lm_summary <- tryCatch(
                summary(lm(clin_var~., data=data_lm)), error=function(e){NULL})
            # Extract regression results
            if (!is.null(lm_summary)) {
                coef <- lm_summary$coefficients[2,]
                tibble::tibble(
                    clin_factor=clin_factor_name, feature=colnames(abundance)[a],
                    method='Linear Regression', beta_coef=coef[1], t_statistic=coef[3],
                    pval=coef[4])
            } else {
                tibble::tibble(
                    clin_factor=clin_factor_name, feature=colnames(abundance)[a],
                    method='Linear Regression', beta_coef=NA, t_statistic=NA, pval=NA)
            }
        })
        return(regression_results)
    }
    lr_table_all <- purrr::map_dfr(
        2:ncol(condition_table), regression_for_factor) %>% # Add adjusted p-values
        dplyr::mutate(
            padj=p.adjust(pval, method=adjust_p_method, n=length(pval)),
            sig_pval=ifelse(pval < p_cutoff, 'yes', 'no'),
            sig_padj=ifelse(padj < p_cutoff, 'yes', 'no') ) %>%
        dplyr::mutate(pval=round(pval, 5), padj=round(padj, 5)) %>% as.data.frame()
    return(lr_table_all=lr_table_all)
}
