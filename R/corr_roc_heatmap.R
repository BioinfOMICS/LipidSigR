#' @title corr_roc_heatmap
#' @description Perform ROC (Receiver Operating Characteristic) analysis for
#' each lipid species or lipid characteristic against specified clinical
#' conditions, and visualize the results as an interactive and static heatmap
#' with hierarchical clustering. Lipids or characteristics with AUC values
#' meeting the specified cutoff are highlighted as significant.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Required when \code{type = 'Char'};
#' must be \code{NULL} when \code{type = 'Sp'}.
#' @param condition_col Character. One or more column names from the group
#' information table representing clinical conditions (e.g., disease status,
#' gene dependency scores) to be used as binary response variables in ROC analysis.
#' @param side_color_char Character. A lipid characteristic used to annotate
#' the side color bar of the heatmap. Must be selected from the common list
#' returned by \code{\link{list_lipid_char}}. Required when \code{type = 'Sp'};
#' must be \code{NULL} when \code{type = 'Char'}.
#' @param auc_cutoff Numeric. Minimum AUC threshold for a feature to be
#' considered significant. Must be a value between 0 and 1. Default is \code{0.5}.
#' @param distfun Character. The correlation method used to compute the distance
#' matrix for hierarchical clustering. Allowed methods include \code{'pearson'},
#' \code{'kendall'}, and \code{'spearman'}. Default is \code{'spearman'}.
#' @param hclustfun Character. The agglomeration method for hierarchical
#' clustering. Must be one of \code{'ward.D'}, \code{'ward.D2'},
#' \code{'single'}, \code{'complete'}, \code{'average'} (= UPGMA),
#' \code{'mcquitty'} (= WPGMA), \code{'median'} (= WPGMC), or
#' \code{'centroid'} (= UPGMC). Default is \code{'average'}.
#' @param transform Character. Method for data transformation applied before
#' analysis. Allowed methods include \code{'none'}, \code{'log10'},
#' \code{'square'}, and \code{'cube'}. Select \code{'none'} to skip
#' transformation. Default is \code{'log10'}.
#' @param type Character. Specifies the unit of correlation analysis:
#' \code{'Sp'} for lipid species-level analysis, or \code{'Char'} for
#' lipid characteristic-level analysis.
#' @return A list containing the following elements:
#' \describe{
#'   \item{all_correlation_result}{A data frame of ROC analysis results for
#'   all features across all clinical conditions, including AUC, 95\% confidence
#'   intervals, and significance flags.}
#'   \item{sig_correlation_result}{A data frame filtered to features whose
#'   AUC meets or exceeds \code{auc_cutoff}.}
#'   \item{interactive_heatmap}{An interactive heatmap widget.}
#'   \item{static_heatmap}{A static heatmap plot.}
#'   \item{heatmap_matrix}{The matrix used for heatmap plotting.}
#' }
#' @export
#' @examples
#' data("corr_data")
#' processed_se <- data_process(
#'     corr_data, exclude_missing=TRUE, exclude_missing_pct=70, replace_na_method='min',
#'     replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
#' result <- corr_roc_heatmap(processed_se, side_color_char='class',
#'     condition_col=c('Sex', 'Smoking'), distfun='spearman',
#'     hclustfun='average', transform='log10', type='Sp')

corr_roc_heatmap <- function(
        processed_se, char=NULL, condition_col=c('Sex', 'Smoking'),
        side_color_char, auc_cutoff=0.5, distfun=c('pearson', 'spearman', 'kendall'),
        hclustfun=c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'),
        transform=c('none', 'log10', 'square', 'cube'), type=c('Sp', 'Char')){

    .check_inputSE(processed_se, metadata_list=NULL)
    ## check for condition_col & adjusted_col
    group_info <- .extract_df(processed_se, type='group')
    .check_df_corrCol(group_info, condition_col, adjusted_col=NULL, cor_type='corr')

    if (!is.numeric(auc_cutoff) | isFALSE(.check_numeric_range(auc_cutoff, 0, 1)) ) {
        stop("auc_cutoff must be a numeric value between 0 and 1.")
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
        abundance_raw <- .extract_df(processed_se, type = "abundance")
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
        abundance_raw <- merge(
            lipid_char_table, char_data, by.x='row.names', by.y='feature') %>%
            dplyr::filter(characteristic == char) %>%
            dplyr::select(-Row.names, -characteristic)
    }
    abundance <- .transform(abundance_raw, transform)
    conCol <- c("sample_name", condition_col)
    condition_table <- .extract_df(processed_se,type='group') %>%
        dplyr::select(all_of(conCol)) %>% as.data.frame()
    uniqueCounts <- sapply(condition_table[-1], function(col) length(unique(col)))
    if(any(uniqueCounts > 2)){
        stop("The 'condition_col' can only select a binary condition.")
    }

    .check_imputation(abundance)

    CHAR <- colnames(abundance)[1]
    colnames(abundance)[1] <- 'feature'
    abundance %<>% tidyr::gather(-feature, key='sample_name', value='value') %>%
        tidyr::spread(key='feature', value='value') %>% dplyr::arrange(sample_name)
    condition_table %<>% dplyr::arrange(sample_name)
    #---------ROC-----------------------------
    feature_names <- colnames(abundance)[-1]
    clin_factors   <- colnames(condition_table)[-1]

    roc_table_all <- lapply(clin_factors, function(clin_factor){
        y <- condition_table[[clin_factor]]
        roc_res <- lapply(feature_names, function(feat){
            ROC_CI <- tryCatch(
                pROC::roc(response=y, predictor=abundance[[feat]], ci=TRUE)$ci,
                error=function(e) NULL)
            if (is.null(ROC_CI)) {
                return(data.frame(feature=feat,
                    roc_auc=NA, roc_auc_lower95=NA, roc_auc_upper95=NA))
            }

            data.frame(feature=feat, roc_auc=ROC_CI[[2]],
                roc_auc_lower95=ROC_CI[[1]], roc_auc_upper95=ROC_CI[[3]])
        })
        dplyr::bind_rows(roc_res) %>%
            dplyr::mutate(clin_factor=clin_factor, method="ROC analysis", .before=1)}) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(sig_auc=ifelse(roc_auc >= auc_cutoff, 'yes','no'))

    roc_table_sig <- roc_table_all %>% dplyr::filter(sig_auc == 'yes')
    heatmap <- suppressWarnings(.corr_heatmap(
        roc_table_all, roc_table_sig, heatmap_col='roc_auc', lipid_char_table, side_color_char,
        significant='roc_auc', distfun, hclustfun, cor_type="lr"))
    if(!is.null(char)){
        colnames(roc_table_all)[2] <- colnames(roc_table_sig)[2] <- char
    }
    return(list(
        all_correlation_result=roc_table_all, sig_correlation_result=roc_table_sig,
        interactive_heatmap=heatmap$in.heatmap, static_heatmap=heatmap$static_hm,
        heatmap_matrix=heatmap$plot.mat))
}
