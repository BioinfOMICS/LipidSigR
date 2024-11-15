#' @title subChar_multiGroup
#' @description This function conducts subgroup lipid characteristic analysis for multi-group data.
#' Lipid species are categorized and summarized into a new lipid abundance table
#' according to two selected lipid characteristics, followed by differential expression analysis.
#' The two chosen characteristics should be either both continuous data or one continuous and one categorical data.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Allowed lipid characteristics include
#' 'Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', and 'FA.OH'.
#' Selected characteristics must be one of 'Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', or 'FA.OH'.
#' @param subChar Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Must be differ from 'char'.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param post_hoc Character. The method to statistic calculation. Allowed method
#' include "One-way ANOVA" and "Kruskal–Wallis test". Default is \code{'One-way ANOVA'}.
#' @param post_hoc_sig Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param post_hoc_p_cutoff Numeric. Significant level. Default is \code{0.05}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("se_multiGroup")
#' processed_se <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' subChar_se <- subChar_multiGroup(
#'     processed_se, char='Total.C', subChar='class', ref_group='ctrl',
#'     post_hoc='One-way ANOVA', post_hoc_sig='pval', post_hoc_p_cutoff=0.05,
#'     transform='log10')


subChar_multiGroup <- function(
        processed_se, char, subChar, ref_group,
        post_hoc=c('One-way ANOVA', 'Kruskal–Wallis test'),
        post_hoc_sig=c('pval', 'padj'), post_hoc_p_cutoff=0.05,
        transform=c('none', 'log10', 'square', 'cube')){

    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    ## List lipid char
    char_list <- list_lipid_char(processed_se)$deChar_list
    if(is.null(char) | isFALSE(char %in% c('Total.C', 'Total.DB', 'Total.OH',
                                           'FA.C', 'FA.DB', 'FA.OH'))){
        stop(paste("The 'char' parameter must be a numeric lipid ",
                   "characteristic, such as 'Total.C'. Please execute",
                   "list_lipid_char function to view the valid lipid characteristics",
                   "you can input."))
    }
    if(char == subChar){
        stop(paste("The 'char' and 'subChar' parameters cannot be the same."))
    }
    if(is.null(subChar) | isFALSE(subChar %in% char_list)){
        stop(paste("The 'subChar' parameter is not in the list of lipid",
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
    if(is.null(transform) | !transform %in% c('none', 'log10', 'square', 'cube')){
        stop(paste("The 'transform' parameter must be one of 'none', 'log10',",
                   "'square', or 'cube'."))
    }

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(processed_se, type="abundance")
    .check_imputation(abundance)
    ## Group information
    group_info <- .extract_df(processed_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'multiple') stop('This function is used for multiple groups. If your data consists of two groups, please use subChar_twoGroup function for analysis.')
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
    lipid_char <- .extract_df(processed_se, type="lipid")

    ## Subset the selected lipid characteristics
    sub_abund <- abundance %>% dplyr::left_join(lipid_char, by='feature') %>%
        dplyr::select('char1'=dplyr::all_of(char), 'char2'=dplyr::all_of(subChar),
                      dplyr::all_of(colnames(abundance))) %>%
        dplyr::filter(!is.na(char1), !is.na(char2)) %>%
        tidyr::separate_rows(char2, sep='\\|') %>%
        tidyr::separate_rows(char1, sep='\\|') %>%
        tidyr::gather(sample_name, abund, -1:-3) %>%
        dplyr::filter_at(dplyr::vars(dplyr::everything()),
                         dplyr::all_vars(!is.na(.))) %>%
        dplyr::group_by(char1, char2, sample_name) %>%
        dplyr::summarise(abund=sum(abund, na.rm=TRUE), .groups='drop') %>%
        tidyr::spread(sample_name, abund) %>%
        dplyr::mutate(feature=paste0(char2, '_', char1)) %>%
        dplyr::select(feature, everything(), -char1, -char2)
    if(nrow(sub_abund) == 0) return(warning(''))
    processed_abund <- sub_abund
    ## Calculate mean and sd of abundance by groups
    mean_abund <- .mean_sd_multiGroup(abundance=sub_abund, group_info)
    ## Data transformation
    sub_abund  <- .transform(abundance=sub_abund, transform)
    ## Add characteristic
    char_abund <- sub_abund %>%
        dplyr::mutate(characteristic=char, .before='feature') %>%
        tidyr::separate(feature, into=c('char2', 'char1'), sep='_', remove=FALSE)
    char2_list <- unique(char_abund$char2)
    ## Two-way ANOVA with interaction effect
    all_test_res <- NULL
    for(i in 1:length(char2_list)){
        CHAR2=ABUND=test_res=NULL
        CHAR2 <- char2_list[i]
        ABUND <- char_abund %>%
            dplyr::filter(char2 == CHAR2) %>%
            dplyr::select(-feature, -char2)
        test_res <- .two_way_anova(abundance=ABUND, group_info, char)
        test_res %<>% dplyr::mutate(subChar=CHAR2, .after='characteristic')
        all_test_res <- data.table::rbindlist(list(all_test_res, test_res))
    }
    ## Post hoc test
    comp_res <- .stat_multiGroup(abundance=sub_abund, group_info, test=post_hoc)
    colnames(comp_res) <- c(
        'feature', 'post_hoc_test', 'post_hoc_statistic', 'post_hoc_pval',
        'post_hoc_negLog10pval', 'post_hoc_padj', 'post_hoc_negLog10padj')
    ## Merge result table
    table_stat <- mean_abund %>%
        tidyr::separate(feature, into=c('sub_feature', 'char1'), sep='_', remove=FALSE) %>%
        dplyr::mutate(sub_characteristic=subChar,
                      characteristic=char, .before='feature') %>%
        dplyr::mutate(method='Two-way ANOVA') %>%
        dplyr::left_join(all_test_res, by=c('characteristic', 'sub_feature'='subChar')) %>%
        dplyr::left_join(comp_res, by='feature') %>%
        dplyr::mutate(
            post_hoc_sig_pval=ifelse(post_hoc_pval < post_hoc_p_cutoff, 'yes', 'no'),
            post_hoc_sig_padj=ifelse(post_hoc_padj < post_hoc_p_cutoff, 'yes', 'no')) %>%
        dplyr::ungroup() %>%
        dplyr::select(-feature) %>%
        dplyr::select(
            sub_characteristic, sub_feature, characteristic, 'feature'='char1',
            dplyr::everything())
    all_deChar_result <- table_stat %>%
        dplyr::select(dplyr::all_of(paste0('post_hoc_sig_', post_hoc_sig)),
                      dplyr::everything()) %>%
        dplyr::arrange(get(paste0('post_hoc_', post_hoc_sig)))
    sig_deChar_result <- all_deChar_result %>%
        dplyr::filter(get(paste0('post_hoc_sig_', post_hoc_sig)) == 'yes')

    ## Transformed abundance
    rownames(sub_abund) <- NULL
    lipid.char <- sub_abund %>% dplyr::select(feature) %>%
        tidyr::separate(feature, into=c('sub_feature', 'char_feature'),
                        sep='_', remove=FALSE) %>%
        dplyr::mutate(sub_characteristic=subChar, characteristic=char) %>%
        dplyr::select(
            feature, sub_characteristic, sub_feature, characteristic, char_feature)
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
                      subChar=subChar,
                      post_hoc_sig=post_hoc_sig,
                      post_hoc_p_cutoff=post_hoc_p_cutoff))
    return(subChar_se=trans_se)
}




