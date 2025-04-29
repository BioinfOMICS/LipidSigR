#' @title deSp_multiGroup
#' @description Compute differentially expressed analysis of multiple groups
#' (independent) to find significant lipid species.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param test Character. The method to use for comparing means. Allowed method include "One-way ANOVA" and
#' "Kruskal-Wallis test". Default is \code{'One-way ANOVA'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. (default: 0.05)
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
#'     normalization='Percentage', transform='log10')
#' deSp_se <- deSp_multiGroup(processed_se, ref_group='ctrl', test='One-way ANOVA',
#'     significant='pval', p_cutoff=0.05, transform='log10')


deSp_multiGroup <- function(processed_se, ref_group,
                            test=c('One-way ANOVA', 'Kruskal-Wallis test'),
                            significant=c('pval', 'padj'), p_cutoff=0.05,
                            transform=c('none', 'log10', 'square', 'cube')){
    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if(is.null(test) | isFALSE(test %in% c('One-way ANOVA', 'Kruskal-Wallis test'))){
        stop("The 'test' parameter must be either 'One-way ANOVA' or 'Kruskal-Wallis test'.")
    }
    if(is.null(significant) | isFALSE(significant %in% c('pval', 'padj'))){
        stop("The 'significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1))){
        stop("The 'p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube'))){
        stop("The 'transform' parameter must be one of 'none', 'log10', 'square' or 'cube'.")
    }
    ## Extract data from SE
    # Lipid abundance
    abundance <- .extract_df(processed_se, type="abundance")
    processed_abund <- abundance
    .check_imputation(abundance)
    if(isTRUE(.check_nor_negative(abundance[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    # Group information
    group_info <- .extract_df(processed_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'multiple') stop('This function is used for multiple groups. If your data consists of two groups, please use deSp_twoGroup function for analysis.')
    group_info %<>% dplyr::mutate(
        original_group_name=group, .before=group,
        group=ifelse(original_group_name == ref_group, 'ctrl',
                     original_group_name))
    if(is.null(ref_group) | isFALSE(ref_group %in% unique(group_info$original_group_name))){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    nSample <- group_info %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(), .groups='drop')
    if(any(nSample$n < 2)){
        stop("Comparing the means of multiple groups requires each group to have at least two samples.")
    }
    lipid_char <- .extract_df(processed_se, type="lipid")
    ## Calculate mean and sd of abuncance by groups
    abund.mean <- .mean_sd_multiGroup(abundance, group_info)
    ## Data transformation
    abundance <- .transform(abundance, transform)
    ## Multi-groups comparison
    comp.res <- .stat_multiGroup(abundance, group_info, test)
    ## Result table
    stat_table <- abund.mean %>% dplyr::left_join(comp.res, by='feature') %>%
        as.data.frame()
    res.tab <- .sig_feature(stat_table, significant, p_cutoff)
    all_deSp_result <- res.tab$all_table %>% dplyr::select(
        dplyr::all_of(paste0('sig_', significant)), dplyr::everything())
    sig_deSp_result <- res.tab$sig_table %>% dplyr::select(
        dplyr::all_of(paste0('sig_', significant)), dplyr::everything())
    if(nrow(sig_deSp_result) == 0){
        warning("This case does not include any significant lipids, which will prevent some subsequent analyses from being performed.")
        sig_deSp_result <- "No significant lipids."
    }
    ## Construct SE objects & Return transformed abundance
    abundance_mat <- abundance %>% dplyr::arrange(feature) %>%
        tibble::column_to_rownames(var="feature")
    lipid_char %<>% dplyr::arrange(feature)
    trans.se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(abundance=as.matrix(abundance_mat)),
        rowData=S4Vectors::DataFrame(lipid_char),
        colData=S4Vectors::DataFrame(group_info),
        metadata=list(all_deSp_result=all_deSp_result,
                      sig_deSp_result=sig_deSp_result,
                      processed_abundance=processed_abund,
                      significant=significant,
                      p_cutoff=p_cutoff,
                      transform=transform))

    return(deSp_se=trans.se)
}





