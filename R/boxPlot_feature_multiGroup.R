#' @title boxPlot_feature_multiGroup
#' @description This function plots an abundance box plot between multiple
#' groups for a user-selected feature (lipid).
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by
#' \code{\link{data_process}}.
#' @param feature Character. A feature (lipid) for plotting.
#' It must be one of the lipids in the lipid abundance data.
#' @param ref_group Character. Group name of the reference group. It must be
#' one of the group names in the group information table's group column.
#' @param test Character. The method for comparing means. Allowed method include
#'  "One-way ANOVA" and "Kruskal-Wallis test". Default is \code{'One-way ANOVA'}.
#' @param post_hoc_sig Character. The p-value to be used for to determine
#' significance, which is only applicable for the post hoc test after the
#' Kruskal-Wallis test. Must be one of "padj" (adjusted p-value) or
#' "pval" (p-value). Default is \code{'pval'}.
#' Note: "One-way ANOVA" can only be "padj", "Kruskal-Wallis test" can be
#' "pval" or "padj"
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", and "cube".
#' Select 'none' to skip data transformation. Default is \code{'log10'}.
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @return Return a list of a box plot and a table.
#' \enumerate{
#' \item static_boxPlot: abundance box plot of the selected feature.
#' \item table_boxplot: table for plotting.
#' \item table_stat: statistical table.
#' }
#' @export
#' @examples
#' data("se_multiGroup")
#' processed_se <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' boxPlot_result <- boxPlot_feature_multiGroup(
#'     processed_se, feature='PE O- 17:0;0_20:3;0', ref_group='ctrl',
#'     test='One-way ANOVA', post_hoc_sig='padj', transform='log10')
boxPlot_feature_multiGroup <- function(
        processed_se, feature, ref_group,
        test=c('One-way ANOVA', 'Kruskal-Wallis test'),
        post_hoc_sig=c('pval', 'padj'),
        transform=c('none', 'log10', 'square', 'cube')){

    ## check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if(is.null(test) |
       isFALSE(test %in% c('One-way ANOVA', 'Kruskal-Wallis test'))){
        stop("The 'test' parameter must be either 'One-way ANOVA' or 'Kruskal-Wallis test'.")
    }
    if(test == 'Kruskal-Wallis test' &
       (is.null(post_hoc_sig) | isFALSE(post_hoc_sig %in% c('pval', 'padj')))){
        stop("The 'post_hoc_sig' parameter must be either 'pval' or 'padj'.")
    }
    if(test == 'One-way ANOVA' & post_hoc_sig == 'pval'){
        post_hoc_sig <- 'padj'
        warning("post_hoc_sig is only for the post hoc test after the ",
                "Kruskal-Wallis test. For the post hoc test after the ",
                "one-way ANOVA, the Tukey's HSD test is used, which only ",
                "returns the adjusted p-value (padj).")
    }
    if(is.null(transform) |
       isFALSE(transform %in% c('none', 'log10', 'square', 'cube'))){
        stop("The transform parameter must be one of 'none', 'log10', 'square' or 'cube'.")
    }

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(processed_se, type="abundance")
    .check_imputation(abundance)
    if(isTRUE(.check_nor_negative(abundance[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    if(!feature %in% abundance$feature){
        stop("'feature' must be one of the lipids in your lipid abundance data; please revise your input.")
    }
    ## Group information
    group_info <- .extract_df(processed_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'multiple'){
        stop('This function is used for multiple groups. If your data consists of two groups, please use boxPlot_feature_twoGroup() for analysis.')
    }
    group_info %<>%
        dplyr::mutate(
            original_group_name=group, .before=group,
            group=ifelse(original_group_name == ref_group, 'ctrl',
                         original_group_name)) %>%
        dplyr::arrange(stringr::str_to_lower(original_group_name)) %>%
        dplyr::arrange(match(group, 'ctrl'))
    if(is.null(ref_group) |
       isFALSE(ref_group %in% unique(group_info$original_group_name))){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    nSample <- group_info %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(), .groups='drop')
    if(any(nSample$n < 2)){
        stop("Comparing the means of multiple groups requires each group to have at least two samples.")
    }
    #### Data transformation ####
    abundance <- .transform(abundance, transform)

    #### Gather the abundance data ####
    sub_abund <- abundance %>% dplyr::filter(feature == .env$feature) %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name')

    #### Multi groups comparison ####
    comp_res <- .box_stat_multiGroup(sub_abund, group_info, ref_group, test)
    #### Box plot ####
    # if(!all(is.na(comp_res$res_table$pval))){
        ## Label of y axis
        if(transform == 'log10'){
            Y <- 'Log10-transformed abundance'
        }else if(transform %in% c('cube', 'sqare')){
            Y <- paste0(stringr::str_to_title(transform),
                        '-transformed abundance')
        }else{
            Y <- 'Processed abundance'
        }
        box.res <- .box_plot_multiGroup(
            sub_abund, group_info, test, post_hoc_sig,
            res_table=comp_res$res_table, post_hoc_table=comp_res$post_hoc_table,
            plot_title=unique(sub_abund$feature), plot_ylabel=Y)
        return(list(static_boxPlot=box.res$box.plot,
                    table_boxplot=box.res$box.tab,
                    table_stat=comp_res$res_table))
}
