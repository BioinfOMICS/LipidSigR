#' @title corr_logR_heatmap
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
#' @param adjusted_col Character. The column names used to extract the adjusted table from
#' the group information table, including additional variables to be incorporated
#' into the algorithm for adjusting confounding effects.
#' @param adjust_p_method Character. The correction method of p-value. Allowed
#' methods include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' and "none". Default is \code{'BH'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{1}.
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
#' @param heatmap_col Character. The value for clustering. Allow method are
#' "ln_OR" and "z_statistic". Default is \code{'ln_OR'}.
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
#'     replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
#' result <- corr_logR_heatmap(
#'     processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
#'     adjusted_col=NULL,
#'     adjust_p_method='BH', significant='pval', p_cutoff=1, side_color_char=NULL,
#'     distfun='spearman', hclustfun='average', heatmap_col='ln_OR',transform='log10', type='Sp')


corr_logR_heatmap <- function(
        processed_se, char=NULL, condition_col, adjusted_col=NULL,
        adjust_p_method='BH', significant='pval', p_cutoff=1, side_color_char,
        distfun='spearman', hclustfun='average', heatmap_col='ln_OR',
        transform='log10', type='Sp'){

    .check_inputSE(processed_se, metadata_list=NULL)
    ## check for condition_col & adjusted_col
    group_info <- .extract_df(processed_se, type='group')
    .check_df_corrCol(group_info, condition_col, adjusted_col, cor_type='lr')


    match.arg(significant, c('pval', 'padj'))
    match.arg(adjust_p_method, c('BH', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BY', 'fdr', 'none'))
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    match.arg(distfun, c('spearman', 'pearson', 'kendall'))
    match.arg(hclustfun, c('average', 'ward.D', 'ward.D2', 'single', 'complete', 'mcquitty', 'median', 'centroid'))
    match.arg(heatmap_col, c('ln_OR', 'z_statistic'))
    match.arg(transform, c('log10', 'square', 'cube', 'none'))
    match.arg(type, c('Sp', 'Char'))

    if (type == "Sp") {
        if (!is.null(char)) {
            stop("If type is set to 'Sp' for running lipid species correlation analysis, char must be set to NULL.")
        }
        if (isFALSE(.check_char(processed_se, side_color_char, type='common'))) {
            stop("Wrong side_color_char input, you can view the available side_color_char list by list_lipid_char function.")
        }
        abundance_raw <- .extract_df(processed_se, type="abundance")
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
        char_data <- .extract_df(char_se, type="abundance")
        lipid_char_table <- .extract_df(char_se, type='lipid')
        abundance_raw <- merge(
            lipid_char_table, char_data, by.x='row.names', by.y='feature') %>%
            dplyr::filter(characteristic == char) %>%
            dplyr::select(-Row.names, -characteristic)
    }
    abundance <- .transform(abundance_raw, transform)
    conCol <- c("sample_name", condition_col)
    condition_table <- .extract_df(processed_se,type='group') %>%
        dplyr::select(all_of(conCol)) %>% as.data.frame()
    .check_imputation(abundance)

    CHAR <- colnames(abundance)[1]
    colnames(abundance)[1] <- 'feature'
    abundance %<>% tidyr::gather(-feature, key='sample_name', value='value') %>%
        tidyr::spread(key='feature', value='value') %>% dplyr::arrange(sample_name)
    condition_table %<>% dplyr::arrange(sample_name)
    if(!is.null(adjusted_col)){
        adjCol <- c("sample_name", adjusted_col)
        adjusted_table <- .extract_df(processed_se,type='group') %>%
            dplyr::select(all_of(adjCol)) %>% as.data.frame() %>%
            dplyr::arrange(sample_name)
    }else{
        adjusted_table <- NULL
    }

    #---------logistic regression-----------------------------
    feature_names <- colnames(abundance)[-1]
    logR_table_all <- lapply(condition_col, function(clin_factor){
        y <- condition_table[[clin_factor]]
        logR_res <- lapply(feature_names, function(feat){
            data_lr <- if (!is.null(adjusted_table)) {
                data.frame(clin_var=y, feature=abundance[[feat]], adjusted_table[-1])
            } else {
                data.frame(clin_var=y,feature=abundance[[feat]])
            }
            lr <- tryCatch(
                stats::glm(formula=clin_var ~ ., family="binomial",
                           data=data_lr, control=list(maxit=50)),
                error=function(e) NULL)
            nullTab <- data.frame(
                feature=feat, ln_OR=NA, OR=NA, OR_lower95=NA, OR_upper95=NA,
                z_statistic=NA, pval=NA)
            if(is.null(lr)){ return(nullTab) }
            coef <- summary(lr)$coefficients

            if(nrow(coef) < 2){ return(nullTab) }

            CI <- tryCatch( confint.default(lr), error=function(e) NULL, warning=function(w) NULL)
            est <- coef[2, 1]

            data.frame(
                feature=feat, ln_OR=est, OR=exp(est),
                OR_lower95=if (!is.null(CI)) exp(CI[2, 1]) else NA,
                OR_upper95=if (!is.null(CI)) exp(CI[2, 2]) else NA,
                z_statistic=coef[2, 3], pval=coef[2, 4])
        })

        dplyr::bind_rows(logR_res) %>%
            dplyr::mutate(clin_factor=clin_factor, method="Logistic Regression",
                padj=p.adjust(pval, method=adjust_p_method))
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            sig_pval=ifelse(pval < p_cutoff, 'yes','no'),
            sig_padj=ifelse(padj < p_cutoff, 'yes','no'))

    logR_table_sig <- logR_table_all %>%
        dplyr::filter(!!rlang::sym(paste0("sig_", significant))=='yes')

    heatmap <- suppressWarnings(.corr_heatmap(
        logR_table_all, logR_table_sig, heatmap_col=heatmap_col, lipid_char_table, side_color_char,
        significant, distfun, hclustfun, cor_type="lr"))
    if(!is.null(char)){
        logR_table_all %<>% dplyr::rename(!!char := feature)
        logR_table_sig %<>% dplyr::rename(!!char := feature)
    }
    return(
        list(
            all_correlation_result=logR_table_all, sig_correlation_result=logR_table_sig,
            interactive_heatmap=heatmap$in.heatmap, static_heatmap=heatmap$static_hm,
            heatmap_matrix=heatmap$plot.mat))
}
