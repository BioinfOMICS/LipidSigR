#' @title deChar_multiGroup
#' @description This function is for multi-group lipid characteristics analysis.
#' It identifies the differentially expressed lipid characteristics based on the user-selected characteristic.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the deChar list
#' returned by \code{\link{list_lipid_char}}.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param post_hoc Character. The method to statistic calculation. Allowed method
#' include "One-way ANOVA" and "Kruskal–Wallis test". Default is \code{'One-way ANOVA'}.
#' @param post_hoc_sig Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param post_hoc_p_cutoff Numeric. Significant level. Default is \code{0.05}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube", "log2". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("se_multiGroup")
#' processed_se <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' char_list <- list_lipid_char(processed_se)$deChar_list
#' print(char_list)
#' deChar_se <- deChar_multiGroup(
#'     processed_se, char='class', ref_group='ctrl', post_hoc='One-way ANOVA',
#'     post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10')


deChar_multiGroup <- function(
        processed_se, char, ref_group,
        post_hoc=c('One-way ANOVA', 'Kruskal–Wallis test'),
        post_hoc_sig=c('pval', 'padj'), post_hoc_p_cutoff=0.05,
        transform=c('none', 'log10', 'square', 'cube', 'log2')){

    #data('lipidChar')
    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    ## List lipid char
    char_list <- list_lipid_char(processed_se)$deChar_list
    if(is.null(char) | isFALSE(char %in% char_list)){
        stop(paste("The 'char' parameter is not in the list of lipid",
                   "characteristics. Please execute list_lipid_char function to view",
                   "the valid lipid characteristics you can input."))
    }
    if(is.null(post_hoc) | isFALSE(post_hoc %in% c('One-way ANOVA', 'Kruskal–Wallis test'))){
        stop("The 'post_hoc' parameter must be either 'One-way ANOVA' or 'Kruskal–Wallis test'.")
    }
    if(is.null(post_hoc_sig) | isFALSE(post_hoc_sig %in% c('pval', 'padj'))){
        stop("The 'post_hoc_sig' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(post_hoc_p_cutoff) | isFALSE(.check_numeric_range(post_hoc_p_cutoff, 0, 1))){
        stop("The 'post_hoc_p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(char %in% lipidChar$characteristic &
       (is.null(transform) | !transform %in% c('none', 'log10', 'square', 'cube'))){
        stop(paste("The 'char' is not a ratio-based characteristic. The",
                   "'transform' parameter must be one of 'none', 'log10',",
                   "'square', or 'cube'."))
    }
    if(!char %in% lipidChar$characteristic &
       (is.null(transform) | !transform %in% c('none', 'log2'))){
        stop(paste("The 'char' is a ratio-based characteristic. The",
                   "'transform' parameter must be either 'none' or 'log2'."))
    }
    if(char %in% lipidChar$characteristic){
        ## convert_sp2char
        char_se <- convert_sp2char(processed_se, transform='none')
    }else{
        ## convert_sp2ratio
        char_se <- convert_sp2ratio(processed_se, transform='none')
    }

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(char_se, type="abundance")
    .check_imputation(abundance)
    ## Group information
    group_info <- .extract_df(char_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'multiple') stop('This function is used for multiple groups. If your data consists of two groups, please use deChar_twoGroup function for analysis.')
    group_info %<>%
        dplyr::mutate(
            original_group_name=group, .before=group,
            group=ifelse(original_group_name == ref_group, 'ctrl',
                         original_group_name)) %>%
        dplyr::arrange(stringr::str_to_lower(original_group_name)) %>%
        dplyr::arrange(match(group, 'ctrl'))
    if(is.null(ref_group) | isFALSE(ref_group %in% unique(group_info$original_group_name))){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    ## Lipid characteristics table
    lipid_char <- .extract_df(char_se, type="lipid")

    ## Subset the selected lipid characteristics
    sub_abund <- merge(lipid_char, abundance, by.x='row.names', by.y='feature') %>%
        dplyr::filter(characteristic == char) %>%
        dplyr::select(-Row.names, -characteristic)
    processed_abund <- sub_abund
    ## Calculate mean and sd of abundance by groups
    mean_abund <- .mean_sd_multiGroup(abundance=sub_abund, group_info)
    ## Data transformation
    sub_abund  <- .transform(abundance=sub_abund, transform)
    ## Add characteristic
    char_abund <- sub_abund %>%
        dplyr::mutate(characteristic=char, .before='feature')

    ## Two-way ANOVA with interaction effect
    test_res <- .two_way_anova(abundance=char_abund, group_info, char=char)
    ## Post hoc test
    comp_res <- .stat_multiGroup(abundance=sub_abund, group_info, test=post_hoc)
    colnames(comp_res) <- c(
        'feature', 'post_hoc_test', 'post_hoc_statistic', 'post_hoc_pval',
        'post_hoc_negLog10pval', 'post_hoc_padj', 'post_hoc_negLog10padj')
    ## Merge result table
    table_stat <- mean_abund %>%
        dplyr::mutate(characteristic=char, .before='feature') %>%
        dplyr::mutate(method='Two-way ANOVA') %>%
        dplyr::left_join(test_res, by='characteristic') %>%
        dplyr::left_join(comp_res, by='feature') %>%
        dplyr::mutate(
            post_hoc_sig_pval=ifelse(post_hoc_pval < post_hoc_p_cutoff, 'yes', 'no'),
            post_hoc_sig_padj=ifelse(post_hoc_padj < post_hoc_p_cutoff, 'yes', 'no'))
    all_deChar_result <- table_stat %>%
        dplyr::select(dplyr::all_of(paste0('post_hoc_sig_', post_hoc_sig)),
                      dplyr::everything()) %>%
        dplyr::arrange(get(paste0('post_hoc_', post_hoc_sig)))
    sig_deChar_result <- all_deChar_result %>%
        dplyr::filter(get(paste0('post_hoc_sig_', post_hoc_sig)) == 'yes')
    if(nrow(sig_deChar_result) == 0){
        warning("This case does not include any significant lipids, which will prevent some subsequent analyses from being performed.")
        sig_deChar_result <- "No significant lipids"
    }
    ## Transformed abundance
    rownames(sub_abund) <- NULL
    lipid.char <- lipid_char %>%
        dplyr::filter(characteristic == char, feature %in% sub_abund$feature) %>%
        dplyr::select(feature, everything())
    rownames(lipid.char) <- NULL
    trans.abund <- sub_abund %>%
        dplyr::arrange(match(feature, lipid.char$feature)) %>%
        tibble::column_to_rownames(var='feature') %>% as.matrix()
    group.info <- group_info %>%
        dplyr::arrange(match(sample_name, colnames(trans.abund)))
    trans_se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(transformed_table=trans.abund),
        colData=S4Vectors::DataFrame(group.info),
        rowData=S4Vectors::DataFrame(lipid.char),
        metadata=list(all_deChar_result=all_deChar_result,
                      sig_deChar_result=sig_deChar_result,
                      processed_abundance=processed_abund,
                      char=char,
                      post_hoc_sig=post_hoc_sig,
                      post_hoc_p_cutoff=post_hoc_p_cutoff))

    return(deChar_se=trans_se)
}
