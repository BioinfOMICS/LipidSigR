#' @title plot_deChar_multiGroup
#' @description This function is for plotting the results of lipid characteristics differential expression analysis.
#' @param deChar_se A SummarizedExperiment object with results computed by \code{\link{deChar_multiGroup}}.
#' @return Return a list with 5 static plots, 4 interactive plots, and 5 data frames.
#' \enumerate{
#' \item static_barPlot: a static bar plot shows the average expression of each sample group and highlights significant differences between groups based on a user-selected characteristic.
#' \item static_barPlot_sqrt: a static bar plot shows the average expression of each sample group and highlights significant differences between groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_linePlot: a static line plot shows the average expression of each sample group and highlights significant differences between the groups based on a user-selected characteristic.
#' \item static_linePlot_sqrt: a static line plot shows the average expression of each sample group and highlights significant differences between the groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_boxPlot: a static box plot of ctrl group and experiment group.
#' \item interactive_barPlot: an interactive bar plot shows the average expression of each sample group and highlights significant differences between groups based on a user-selected characteristic.
#' \item interactive_barPlot_sqrt: an interactive bar plot shows the average expression of each sample group and highlights significant differences between groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item interactive_linePlot: an interactive line plot shows the average expression of each sample group and highlights significant differences between the groups based on a user-selected characteristic.
#' \item interactive_linePlot_sqrt: an interactive line plot shows the average expression of each sample group and highlights significant differences between the groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item table_barPlot: table for plotting bar plots
#' \item table_linePlot: table for plotting line plots
#' \item table_boxPlot: table for plotting box plots
#' \item table_char_index: table with the value calculated by the weighted average of lipid characteristics abundance
#' \item table_index_stat: table with statistics of control and experiment groups
#' }
#' @export
#' @examples
#' data("se_multiGroup")
#' processed_se <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deChar_se <- deChar_multiGroup(
#'     processed_se, char='Total.C', ref_group='ctrl', post_hoc='One-way ANOVA',
#'     post_hoc_sig='pval', post_hoc_p_cutoff=0.05, transform='log10')
#' deChar_plot <- plot_deChar_multiGroup(deChar_se)

plot_deChar_multiGroup <- function(deChar_se){

    #### Check parameter ####
    .check_de_outputSE(deChar_se, de_type="deChar")

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- S4Vectors::metadata(deChar_se)$processed_abundance
    .check_imputation(abundance)
    ## Group information
    group_info <- .extract_df(deChar_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'multiple'){
        stop(paste("This function is used for multiple groups. If your data",
                   "consists of two groups, please use plot_deChar_twoGroup",
                   "for analysis."))
    }
    ref_group <- unique(group_info$original_group_name[which(group_info$group == 'ctrl')])
    group_info %<>%
        dplyr::arrange(stringr::str_to_lower(original_group_name)) %>%
        dplyr::arrange(match(group, 'ctrl'))
    ## Meta data
    all_deChar_result <- S4Vectors::metadata(deChar_se)$all_deChar_result
    post_hoc <- unique(all_deChar_result$post_hoc_test)
    post_hoc_sig <- S4Vectors::metadata(deChar_se)$post_hoc_sig
    char <- S4Vectors::metadata(deChar_se)$char

    #### Create bar plot & line plot table ####
    plot_tab <- .deChar_plot_tab_multiGroup(
        stat_tab=all_deChar_result, group_info, post_hoc_sig)
    ## Bar plot
    bar_res <- .deChar_bar_multiGroup(plot_tab, char, star_position=5)
    ## Continuous characteristics
    if(char %in% c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH')){
        ## Line plot
        line_res <- .deChar_line_multiGroup(plot_tab, char, star_position=5)
        ## Combined char into one
        index_abund <- .char_index(abundance, char)
        ## Calculate mean and sd of abundance by groups
        index_mean <- .mean_sd_multiGroup(abundance=index_abund, group_info)
        ## Statistical test
        index_stat <- .stat_multiGroup(abundance=index_abund, group_info,
                                       test=post_hoc)
        index_stat_tab <- index_mean %>% dplyr::left_join(index_stat, by='feature')
        #### Box plot ####
        ## Gather the abundance data
        sub_abund <- index_abund %>%
            tidyr::gather(sample_name, abund, -1) %>%
            dplyr::left_join(group_info, by='sample_name')
        ## Multi groups comparison
        comp_res <- .box_stat_multiGroup(sub_abund, group_info, ref_group, test=post_hoc)
        ## Box plot
        if(!all(is.na(comp_res$res_table$pval))){
            signif <- ifelse(post_hoc == 'One-way ANOVA', 'padj', 'pval')
            box_res <- .box_plot_multiGroup(
                sub_abund, group_info, test=post_hoc, post_hoc_sig=signif,
                res_table=comp_res$res_table, post_hoc_table=comp_res$post_hoc_table,
                plot_title=NULL, plot_ylabel=paste0(char, ' index'))
        }else{
            box_res <- NULL
        }
        return(list(static_barPlot=bar_res$bar.plot,
                    static_barPlot_sqrt=bar_res$bar.sqrt.plot,
                    static_linePlot=line_res$line.plot,
                    static_linePlot_sqrt=line_res$line.sqrt.plot,
                    static_boxPlot=box_res$box.plot,
                    interactive_barPlot=bar_res$bar.plotly,
                    interactive_barPlot_sqrt=bar_res$bar.sqrt.plotly,
                    interactive_linePlot=line_res$line.plotly,
                    interactive_linePlot_sqrt=line_res$line.sqrt.plotly,
                    table_barPlot=bar_res$bar.tab,
                    table_linePlot=line_res$line.tab,
                    table_boxPlot=box_res$box.tab,
                    table_char_index=index_abund,
                    table_index_stat=index_stat_tab))
    }else{
        message(paste("The 'char' parameter is not a numeric characteristic,",
                      "so it is unable to perform line plot and box plot."))
        return(list(static_barPlot=bar_res$bar.plot,
                    static_barPlot_sqrt=bar_res$bar.sqrt.plot,
                    interactive_barPlot=bar_res$bar.plotly,
                    interactive_barPlot_sqrt=bar_res$bar.sqrt.plotly,
                    table_barPlot=bar_res$bar.tab))
    }
}






