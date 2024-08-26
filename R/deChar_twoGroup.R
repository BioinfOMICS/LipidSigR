#' @title deChar_twoGroup
#' @description This function is for two-group lipid characteristics analysis.
#' It identifies the differentially expressed lipid characteristics based on the user-selected characteristic.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the deChar list
#' returned by \code{\link{list_lipid_char}}.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param test Character. The method to use for comparing means.
#' Allowed method include "t-test" and "Wilcoxon test". Default is \code{'t-test'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{0.05}.
#' @param FC_cutoff Numeric. Significance of the fold-change. Default is \code{1}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube", "log2". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' char_list <- list_lipid_char(processed_se)$deChar_list
#' names(char_list) <- NULL
#' deChar_se <- deChar_twoGroup(processed_se, char=char_list[14], ref_group="ctrl",
#'     test='t-test', significant="pval", p_cutoff=0.05,
#'     FC_cutoff=1, transform='log10')
#'
deChar_twoGroup <- function(
        processed_se, char, ref_group, test=c('t-test', 'Wilcoxon test'),
        significant=c("pval", "padj"),
        p_cutoff=0.05, FC_cutoff=1, transform=c('none', 'log10', 'square', 'cube', 'log2')){

    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    ## List lipid char
    char_list <- list_lipid_char(processed_se)$deChar_list
    if(is.null(char) | isFALSE(char %in% char_list)){
        stop(paste("The 'char' parameter is not in the list of lipid",
                   "characteristics. Please execute list_lipid_char() to view",
                   "the valid lipid characteristics you can input."))
    }
    if (is.null(test) | isFALSE(test %in% c('t-test', 'Wilcoxon test')) ) {
        stop("test must be one of 't-test' or 'Wilcoxon test''.")
    }
    if (is.null(significant) | isFALSE(significant %in% c('pval', 'padj')) ) {
        stop("significant must be one of 'pval' or 'padj'.")
    }
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    if (!is.numeric(FC_cutoff) | isFALSE(.check_numeric_range(FC_cutoff, 1, 8)) ) {
        stop("FC_cutoff must be a numeric value between 1 and 8.")
    }
    if(char %in% lipidChar$characteristic &
       (is.null(transform) | !transform %in% c('none', 'log10', 'square', 'cube'))){
        stop("The 'char' is not a ratio-based characteristic. The 'transform'
             parameter must be one of 'none', 'log10', 'square', or 'cube'.")
    }
    if(!char %in% lipidChar$characteristic &
       (is.null(transform) | !transform %in% c('none', 'log2'))){
        stop("The 'char' is a ratio-based characteristic. The 'transform'
             parameter must be either 'none' or 'log2'.")
    }
    if(char %in% lipidChar$characteristic){
        ## convert_sp2char
        char_se <- convert_sp2char(processed_se, transform='none')
    }else{
        ## convert_sp2ratio
        char_se <- convert_sp2ratio(processed_se, transform='none')
    }

    #### Extract data from SE ####
    abundance <- .extract_df(char_se, type="abundance")
    group_info_raw <- .extract_df(char_se, type="group")
    .check_imputation(abundance)
    if (isTRUE(.check_nor_negative(abundance[-1])) ) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    if (is.null(ref_group) | isFALSE(ref_group %in% unique(group_info_raw$group)) ){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    ## Group information
    if (.check_nGroup(group_info_raw) != 'two') {
        stop('This function is used for two groups. If your data consists of multi groups, please use deChar_nultiGroup() for analysis.')
    }
    ## add original group names
    group_info <- group_info_raw %>%
        dplyr::mutate("original_group_name"=group_info_raw$group, .before=group)
    group_info[group_info$original_group_name==ref_group, "group"] <- "ctrl"
    group_info[group_info$original_group_name!=ref_group, "group"] <- "exp"
    ## paired sample
    paired_sample <- ifelse(isFALSE(all(is.na(group_info$pair))), TRUE, FALSE)
    ## mean_ctrl, mean_exp, method, FC, log2FC
    abundance_tab <- .mean_sd_twoGroup(abundance, group_info)

    ## Lipid characteristics table
    lipid_char <- .extract_df(char_se, type="lipid")

    ## Subset the selected lipid characteristics
    sub_abund <- merge(lipid_char, abundance, by.x='row.names', by.y='feature') %>%
        dplyr::filter(characteristic == char) %>%
        dplyr::select(-Row.names, -characteristic)
    processed_abund <- sub_abund
    ## Calculate mean and sd of abundance by groups
    mean_abund <- .mean_sd_twoGroup(abundance=sub_abund, group_info)
    ## Data transformation
    sub_abund  <- .transform(abundance=sub_abund, transform)
    ## Add characteristic
    char_abund <- sub_abund %>%
        dplyr::mutate(characteristic=char, .before='feature')
    ## Two-way ANOVA with interaction effect
    test_res <- .two_way_anova(abundance=char_abund, group_info, char=char)
    ## Post hoc test
    comp_res <- .stat_twoGroup(abundance=sub_abund, group_info, test=test, paired_sample=paired_sample)
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
            post_hoc_sig_pval=ifelse(post_hoc_pval < p_cutoff, 'yes', 'no'),
            post_hoc_sig_padj=ifelse(post_hoc_padj < p_cutoff, 'yes', 'no'))
    all_deChar_result <- table_stat %>%
        dplyr::select(dplyr::all_of(paste0('post_hoc_sig_', significant)),
                      dplyr::everything()) %>%
        dplyr::arrange(get(paste0('post_hoc_', significant)))
    sig_deChar_result <- all_deChar_result %>%
        dplyr::filter(get(paste0('post_hoc_sig_', significant)) == 'yes')
    ## Transformed abundance
    rownames(sub_abund) <- NULL
    lipid.char <- lipid_char %>%
        dplyr::filter(characteristic==char, feature %in% sub_abund$feature) %>%
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
                      significant=significant,
                      p_cutoff=p_cutoff,
                      FC_cutoff=FC_cutoff))
    return(deChar_se=trans_se)
}

#' @title plot_deChar_twoGroup
#' @description This function is for plotting the results of lipid characteristics differential expression analysis.
#' @param deChar_se A SummarizedExperiment object with results computed by \code{\link{deChar_twoGroup}}.
#' @return Return a list with 5 static plots, 5 interactive plots, and 5 data frames.
#' \enumerate{
#' \item static_barPlot: a static bar plot shows the average expression of each sample group and highlights significant differences between two groups based on a user-selected characteristic.
#' \item static_barPlot_sqrt: a static bar plot shows the average expression of each sample group and highlights significant differences between two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_linePlot: a static line plot shows the average expression of each sample group and highlights significant differences between the two groups based on a user-selected characteristic.
#' \item static_linePlot_sqrt: a static line plot shows the average expression of each sample group and highlights significant differences between the two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_boxPlot: a static box plot of ctrl group and experiment group.
#' \item interactive_barPlot: an interactive bar plot shows the average expression of each sample group and highlights significant differences between two groups based on a user-selected characteristic.
#' \item interactive_barPlot_sqrt: an interactive bar plot shows the average expression of each sample group and highlights significant differences between two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item interactive_linePlot: an interactive line plot shows the average expression of each sample group and highlights significant differences between the two groups based on a user-selected characteristic.
#' \item interactive_linePlot_sqrt: an interactive line plot shows the average expression of each sample group and highlights significant differences between the two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item interactive_boxPlot: an interactive box plot of ctrl group and experiment group.
#' \item table_barPlot: table for plotting bar plots
#' \item table_linePlot: table for plotting line plots
#' \item table_boxPlot: table for plotting box plots
#' \item table_char_index: table with the value calculated by the weighted average of lipid characteristics abundance
#' \item table_index_stat: table with statistics of control and experiment groups
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' deChar_se <- deChar_twoGroup(processed_se, char="Total.C", ref_group="ctrl",
#'     test='t-test', significant="padj", p_cutoff=0.05,
#'     FC_cutoff=1, transform='log10')
#' res_plot <- plot_deChar_twoGroup(deChar_se)
plot_deChar_twoGroup<- function(deChar_se){
    #### Check parameter ####
    .check_inputSE(
        deChar_se, metadata_list=c(
            'all_deChar_result', 'sig_deChar_result', 'processed_abundance',
            'char', 'significant', 'p_cutoff', 'FC_cutoff'))
    abundance <- S4Vectors::metadata(deChar_se)$processed_abundance
    .check_imputation(abundance)
    lipid_char_table <- .extract_df(deChar_se, type="lipid")
    ## Group information
    group_info <- .extract_df(deChar_se, type="group")
    if(.check_nGroup(group_info) != 'two'){
        stop(paste("This function is used for two groups. If your data",
                   "consists of two groups, please use plot_deChar_multiGroup()",
                   "for analysis."))
    }
    ref_group <- unique(group_info$original_group_name[which(group_info$group == 'ctrl')])
    paired_sample <- ifelse(isFALSE(all(is.na(group_info$pair))), TRUE, FALSE)
    group_info %<>%
        dplyr::arrange(stringr::str_to_lower(original_group_name)) %>%
        dplyr::arrange(match(group, 'ctrl'))
    ## Meta data
    all_deChar_result <- S4Vectors::metadata(deChar_se)$all_deChar_result
    post_hoc_test <- unique(all_deChar_result$post_hoc_test)
    significant <- S4Vectors::metadata(deChar_se)$significant
    char <- S4Vectors::metadata(deChar_se)$char
    p_cutoff <- S4Vectors::metadata(deChar_se)$p_cutoff
    FC_cutoff <- S4Vectors::metadata(deChar_se)$FC_cutoff
    ## plot table
    plot_tab <- .table_deChar_twoGroup(stat_tab=all_deChar_result, group_info, significant)
    ## Bar plot
    bar_plot <- .barPlot_deChar_twoGroup(plot_tab, char, split_class=NULL)

    if(char %in% c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH')){
        ## Line plot
        line_plot <- .lintPlot_deChar_twoGroup(plot_tab, char, split_class=NULL)
        ## Combined char into one
        index_abund <- .char_index(abundance, char)
        ## Calculate mean and sd of abundance by groups
        index_mean <- .mean_sd_twoGroup(abundance=index_abund, group_info)
        ## Statistical test
        index_stat <- .stat_twoGroup(abundance=index_abund, group_info, paired_sample, test=post_hoc_test)
        index_stat_tab <- index_mean %>% dplyr::left_join(index_stat, by='feature')
        Combine_char_result_table <- index_stat_tab
        Combine_char_result_table <- Combine_char_result_table %>% dplyr::mutate(sig=ifelse((significant < p_cutoff) & (abs(log2FC) > log2(FC_cutoff)),'yes', 'no'))
        Combine_char_result_table <- Combine_char_result_table %>% dplyr::mutate(sig=ifelse(is.na(sig), 'no', sig ))
        #### Box plot ####
        ## Gather the abundance data
        sub_abund <- index_abund %>%
            tidyr::gather(sample_name, abund, -1) %>%
            dplyr::left_join(group_info, by='sample_name')
        box_plot <- .boxPlot_deChar_twoGroup(Combine_char_result_table, boxTab=sub_abund, char, split_class=NULL)

        return(list(static_barPlot=bar_plot$barPlot,
                    static_barPlot_sqrt=bar_plot$barPlot_sqrt,
                    static_linePlot=line_plot$linePlot,
                    static_linePlot_sqrt=line_plot$linePlot_sqrt,
                    static_boxPlot=box_plot$boxPlot,
                    interactive_barPlot=bar_plot$in_barPlot,
                    interactive_barPlot_sqrt=bar_plot$in_barPlot_sqrt,
                    interactive_linePlot=line_plot$in_linePlot,
                    interactive_linePlot_sqrt=line_plot$in_linePlot_sqrt,
                    interactive_boxPlot=box_plot$in_boxPlot,
                    table_barPlot=plot_tab,
                    table_linePlot=plot_tab,
                    table_boxPlot=sub_abund,
                    table_char_index=index_abund,
                    table_index_stat=index_stat_tab))
    }else{
        message("The 'char' parameter is not a numeric characteristic so it is unable to perform line plot and box plot.")
        return(list(static_barPlot=bar_plot$barPlot,
                    static_barPlot_sqrt=bar_plot$barPlot_sqrt,
                    interactive_barPlot=bar_plot$in_barPlot,
                    interactive_barPlot_sqrt=bar_plot$in_barPlot_sqrt,
                    table_barPlot=plot_tab))
    }
}
