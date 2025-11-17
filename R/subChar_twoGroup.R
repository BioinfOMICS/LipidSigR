#' @title subChar_twoGroup
#' @description This function conducts subgroup lipid characteristic analysis for two-group data.
#' Lipid species are categorized and summarized into a new lipid abundance table
#' according to two selected lipid characteristics, followed by differential expression analysis.
#' The two chosen characteristics should be either both continuous data or one continuous and one categorical data.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Allowed lipid characteristics include
#' 'Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', and 'FA.OH'.
#' @param subChar Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Must be differ from 'char'.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param test A character string indicating which method to be used for
#' comparing means. Allowed method include \bold{"t-test"} and
#' \bold{"Wilcoxon test"}. Default is \code{'t-test'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{0.05}.
#' @param FC_cutoff Numeric. Significance of the fold-change. Default is \code{1}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' subChar_se <- subChar_twoGroup(processed_se, char="Total.C", subChar="class",
#'     ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
#'     FC_cutoff=1, transform='log10')

subChar_twoGroup <- function(
    processed_se, char, subChar, ref_group, test=c('t-test', 'Wilcoxon test'),
    significant=c("pval", "padj"), p_cutoff=0.05, FC_cutoff=1,
    transform=c('none', 'log10', 'square', 'cube')) {

    #data('lipidChar')
    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    ## List lipid char
    char_list <- list_lipid_char(processed_se)$deChar_list
    if(is.null(char) | isFALSE(char %in% c('Total.C', 'Total.DB', 'Total.OH',
                                           'FA.C', 'FA.DB', 'FA.OH'))){
        stop("The 'char' parameter must be a numeric lipid characteristic, such as 'Total.C'. Please execute list_lipid_char function to view the valid lipid characteristics you can input.")
    }
    if(char == subChar){
        stop("The 'char' and 'subChar' parameters cannot be the same.")
    }
    if(is.null(subChar) | isFALSE(subChar %in% char_list)){
        stop("The 'subChar' parameter is not in the list of lipid characteristics. Please execute list_lipid_char function to view the valid lipid characteristics you can input.")
    }
    if(is.null(test) | isFALSE(test %in% c('t-test', 'Wilcoxon test'))){
        stop("The 'test' parameter must be either 't-test', or 'Wilcoxon test'.")
    }
    if(is.null(significant) | isFALSE(significant %in% c('pval', 'padj'))){
        stop("The 'significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1))){
        stop("The 'p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(is.null(transform) | !transform %in% c('none', 'log10', 'square', 'cube')){
        stop("The 'transform' parameter must be one of 'none', 'log10', 'square', or 'cube'.")
    }

    #### Extract data from SE ####
    abundance <- .extract_df(processed_se, type="abundance")
    group_info_raw <- .extract_df(processed_se, type="group")
    .check_imputation(abundance)
    if (isTRUE(.check_nor_negative(abundance[-1])) ) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    if (is.null(ref_group) | isFALSE(ref_group %in% unique(group_info_raw$group)) ){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    ## Group information
    if (.check_nGroup(group_info_raw) != 'two') {
        stop('This function is used for two groups. If your data consists of multi groups, please use deChar_multiGroup functionfor analysis.')
    }
    ## add original group names
    group_info <- group_info_raw %>%
        dplyr::mutate("original_group_name"=group_info_raw$group, .before=group)
    group_info[group_info$original_group_name==ref_group, "group"] <- "ctrl"
    group_info[group_info$original_group_name!=ref_group, "group"] <- "exp"
    ## paired sample
    paired_sample <- ifelse(isFALSE(all(is.na(group_info$pair))), TRUE, FALSE)
    ## Lipid characteristics table
    lipid_char <- .extract_df(processed_se, type="lipid")
    ## Subset the selected lipid characteristics
    sub_abund <- abundance %>% dplyr::left_join(lipid_char, by='feature') %>%
        dplyr::select('char1'=dplyr::all_of(char), 'char2'=dplyr::all_of(subChar),
                      dplyr::all_of(colnames(abundance))) %>%
        dplyr::filter(!is.na(char1), !is.na(char2)) %>%
        tidyr::separate_rows(char2, sep='\\|') %>%
        tidyr::separate_rows(char1, sep='\\|') %>%
        tidyr::gather(sample_name, abund, -seq_len(3)) %>%
        dplyr::filter_at(dplyr::vars(dplyr::everything()),
                         dplyr::all_vars(!is.na(.))) %>%
        dplyr::group_by(char1, char2, sample_name) %>%
        dplyr::summarise(abund=sum(abund, na.rm=TRUE), .groups='drop') %>%
        tidyr::spread(sample_name, abund) %>%
        dplyr::mutate(feature=paste0(char2, '_', char1)) %>%
        dplyr::select(feature, everything(), -char1, -char2)
    if(nrow(sub_abund) == 0) {
        stop("Insufficient data for further analysis, please re-select the char and subChar.")
    }
    processed_abund <- sub_abund
    ## Calculate mean and sd of abundance by groups
    mean_abund <- .mean_sd_twoGroup(abundance=sub_abund, group_info)
    ## Data transformation
    sub_abund  <- .transform(abundance=sub_abund, transform)
    ## Add characteristic
    char_abund <- sub_abund %>%
        dplyr::mutate(characteristic=char, .before='feature') %>%
        tidyr::separate(feature, into=c('char2', 'char1'), sep='_', remove=FALSE)
    char2_list <- unique(char_abund$char2)
    ## Two-way ANOVA with interaction effect
    all_test_res <- NULL
    for(i in seq_len(length(char2_list))){
        CHAR2 <- ABUND <- test_res <- NULL
        CHAR2 <- char2_list[i]
        ABUND <- char_abund %>%
            dplyr::filter(char2 == CHAR2) %>%
            dplyr::select(-feature, -char2)
        test_res <- .two_way_anova(abundance=ABUND, group_info, char)
        test_res %<>% dplyr::mutate(subChar=CHAR2, .after='characteristic')
        all_test_res <- data.table::rbindlist(list(all_test_res, test_res))
    }
    ## Post hoc test
    comp_res <- .stat_twoGroup(abundance=sub_abund, group_info, paired_sample, test)
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
            post_hoc_sig_pval=ifelse(post_hoc_pval < p_cutoff, 'yes', 'no'),
            post_hoc_sig_padj=ifelse(post_hoc_padj < p_cutoff, 'yes', 'no')) %>%
        dplyr::ungroup() %>%
        dplyr::select(-feature) %>%
        dplyr::select(
            sub_characteristic, sub_feature, characteristic, 'feature'='char1',
            dplyr::everything())
    all_deChar_result <- table_stat %>%
        dplyr::select(dplyr::all_of(paste0('post_hoc_sig_', significant)),
                      dplyr::everything()) %>%
        dplyr::arrange(get(paste0('post_hoc_', significant)))
    sig_deChar_result <- all_deChar_result %>%
        dplyr::filter(get(paste0('post_hoc_sig_', significant)) == 'yes')

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
                      significant=significant,
                      p_cutoff=p_cutoff,
                      FC_cutoff=FC_cutoff))
    return(subChar_se=trans_se)
}

#' @title plot_subChar_twoGroup
#' @description This function is for plotting the results of subgroup lipid characteristics differential expression analysis.
#' @param subChar_se A SummarizedExperiment object with results computed by \code{\link{subChar_twoGroup}}.
#' @param subChar_feature Character. A feature selected by users from \code{subChar} to visualize the specific plot for the selected category of that characteristic.
#' For example, if subChar is 'class' and subChar_feature is 'Cer', the resulting plots will display data for 'Cer' within the 'class' category.
#' @return Return a list with 5 static plots, 5 interactive plots, and 5 data frames.
#' \enumerate{
#' \item static_barPlot: a static bar plot shows the average expression of each sample group and highlights significant
#' differences between two groups based on a user-selected characteristic.
#' \item static_barPlot_sqrt: a static bar plot shows the average expression of each sample
#' group and highlights significant differences between two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_linePlot: a static line plot shows the average expression of each sample group and
#' highlights significant differences between the two groups based on a user-selected characteristic.
#' \item static_linePlot_sqrt: a static line plot shows the average expression of each sample group and highlights
#' significant differences between the two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item static_boxPlot: a static box plot of ctrl group and experiment group.
#' \item interactive_barPlot: an interactive bar plot shows the average expression of each sample group
#' and highlights significant differences between two groups based on a user-selected characteristic.
#' \item interactive_barPlot_sqrt: an interactive bar plot shows the average expression of each sample
#' group and highlights significant differences between two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
#' \item interactive_linePlot: an interactive line plot shows the average expression of each
#' sample group and highlights significant differences between the two groups based on a user-selected characteristic.
#' \item interactive_linePlot_sqrt: an interactive line plot shows the average expression of
#' each sample group and highlights significant differences between the two groups based on a user-selected characteristic. \emph{NOTE: the y axis is sqrt-scaled}.
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
#'     normalization='Percentage', transform='log10')
#' subChar_se <- subChar_twoGroup(processed_se, char="Total.C", subChar="class",
#'     ref_group="ctrl", test='t-test', significant="pval", p_cutoff=0.05,
#'     FC_cutoff=1, transform='log10')
#' subChar_feature_list <- unique(
#'     extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)
#' res_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")
plot_subChar_twoGroup <- function(subChar_se, subChar_feature){
    #### Check parameter ####
    .check_inputSE(
        subChar_se, metadata_list=c(
            'all_deChar_result', 'sig_deChar_result', 'processed_abundance',
            'char', 'subChar', 'significant', 'p_cutoff', 'FC_cutoff'))
    abundance <- S4Vectors::metadata(subChar_se)$processed_abundance
    .check_imputation(abundance)
    lipid_char_table <- .extract_df(subChar_se, type="lipid")
    ## Group information
    group_info <- .extract_df(subChar_se, type="group")
    if(.check_nGroup(group_info) != 'two'){
        stop("This function is used for two groups. If your data", "consists of two groups, please use plot_deChar_multiGroup function", "for analysis.")
    }
    ref_group <- unique(group_info$original_group_name[which(group_info$group == 'ctrl')])
    paired_sample <- ifelse(isFALSE(all(is.na(group_info$pair))), TRUE, FALSE)
    group_info %<>%
        dplyr::arrange(stringr::str_to_lower(original_group_name)) %>%
        dplyr::arrange(match(group, 'ctrl'))
    ## Meta data
    all_deChar_result <- S4Vectors::metadata(subChar_se)$all_deChar_result
    if(is.null(subChar_feature) | isFALSE(subChar_feature %in% unique(all_deChar_result$sub_feature))){
        stop("The 'subChar_feature' parameter must be one of the valid subChar features.")
    }
    post_hoc_test <- unique(all_deChar_result$post_hoc_test)
    significant <- S4Vectors::metadata(subChar_se)$significant
    char <- S4Vectors::metadata(subChar_se)$char
    subChar <- S4Vectors::metadata(subChar_se)$subChar
    p_cutoff <- S4Vectors::metadata(subChar_se)$p_cutoff
    FC_cutoff <- S4Vectors::metadata(subChar_se)$FC_cutoff

    #### Subset the data ####
    subChar_abund <- abundance %>%
        dplyr::filter(stringr::str_starts(feature, paste0(subChar_feature, '_'))) %>%
        dplyr::mutate(feature=stringr::str_remove(feature, paste0(subChar_feature, '_')))
    subChar_res <- all_deChar_result %>%
        dplyr::filter(sub_feature == subChar_feature)
    #### Create bar plot & line plot table ####
    plot_tab <- .table_deChar_twoGroup(stat_tab=subChar_res, group_info, significant)
    ## Bar plot
    bar_res <- .barPlot_deChar_twoGroup(plot_tab, char, split_class=subChar_feature)
    ## Line plot
    line_res <- .lintPlot_deChar_twoGroup(plot_tab, char, split_class=subChar_feature)
    ## Combined char into one
    index_abund <- .char_index(abundance=subChar_abund, char)
    ## Calculate mean and sd of abundance by groups
    index_mean <- .mean_sd_twoGroup(abundance=index_abund, group_info)
    ## Statistical test
    index_stat <- .stat_twoGroup(abundance=index_abund, group_info, paired_sample, test=post_hoc_test)
    index_stat_tab <- index_mean %>% dplyr::left_join(index_stat, by='feature')
    Combine_char_result_table <- index_stat_tab
    Combine_char_result_table <- Combine_char_result_table %>% dplyr::mutate(sig=ifelse((significant < p_cutoff) & (abs(log2FC) > log2(FC_cutoff)),'yes', 'no'))
    Combine_char_result_table <- Combine_char_result_table %>% dplyr::mutate(sig=ifelse(is.na(sig), 'no', sig ))
    ## Gather the abundance data
    sub_abund <- index_abund %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name')
    box_res <- .boxPlot_deChar_twoGroup(Combine_char_result_table, boxTab=sub_abund, char, split_class=subChar_feature)

    return(list(static_barPlot=bar_res$barPlot,
                static_barPlot_sqrt=bar_res$barPlot_sqrt,
                static_linePlot=line_res$linePlot,
                static_linePlot_sqrt=line_res$linePlot_sqrt,
                static_boxPlot=box_res$boxPlot,
                interactive_barPlot=bar_res$in_barPlot,
                interactive_barPlot_sqrt=bar_res$in_barPlot_sqrt,
                interactive_linePlot=line_res$in_linePlot,
                interactive_linePlot_sqrt=line_res$in_linePlot_sqrt,
                interactive_boxPlot=box_res$in_boxPlot,
                table_barPlot=plot_tab,
                table_linePlot=plot_tab,
                table_boxPlot=sub_abund,
                table_char_index=index_abund,
                table_index_stat=index_stat_tab))
}
