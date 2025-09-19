#' @title boxPlot_feature_twoGroup
#' @description This function plots an abundance box plot between two groups
#' for a user-selected feature (lipid).
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by
#' \code{\link{data_process}}.
#' @param feature Character. A feature (lipid) for plotting. It must be one of
#' the lipids in the lipid abundance data.
#' @param ref_group Character. Group name of the reference group. It must be
#' one of the group names in the group information table's group column.
#' @param test Character. The method for comparing means. Allowed method include
#' "t-test" and "Wilcoxon test". Default is \code{'t-test'}
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", and "cube". Select 'none' to skip data
#' transformation. Default is \code{'log10'}.
#' @return Return a list of a box plot and a table.
#' \enumerate{
#' \item static_boxPlot: abundance box plot of the selected feature.
#' \item table_boxplot: table for plotting.
#' \item table_stat: statistical table.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' boxPlot_result <- boxPlot_feature_twoGroup(
#'     processed_se, feature='Cer 38:1;2', ref_group='ctrl', test='t-test',
#'     transform='log10')
boxPlot_feature_twoGroup <- function(
        processed_se, feature, ref_group, test=c('t-test', 'Wilcoxon test'),
        transform=c('none', 'log10', 'square', 'cube')){

    ## check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if (is.null(test) | isFALSE(test %in% c('t-test', 'Wilcoxon test')) ) {
        stop("The 'test' parameter must be either 't-test' or 'Wilcoxon test'.")
    }
    if (is.null(transform) |
        isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
        stop("The 'transform' parameter must be one of 'none', 'log10',",
             " 'square' or 'cube'.")
    }
    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(processed_se, type="abundance")
    .check_imputation(abundance)
    if(isTRUE(.check_nor_negative(abundance[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended
             choosing a different normalization method during data processing.")
    }
    if(!feature %in% abundance$feature){
        stop("'feature' must be one of the lipids in your lipid abundance data; ",
        "please revise your input.")
    }
    ## Group information
    group_info <- .extract_df(processed_se, type="group")
    group_type <- .check_nGroup(group_info)
    if(group_type != 'two'){
        stop('This function is used for two groups. If your data consists of ',
             'multiple groups, please use "boxPlot_feature_multiGroup"',
             ' for analysis.')
    }
    if(sum(group_info$group == 'ctrl') < 2 |
       sum(group_info$group == 'exp') < 2 ){
        stop("Comparing the means of two independent groups requires each group
             to have at least two samples.")
    }
    group_info %<>% dplyr::mutate(
        original_group_name=group, .before=group,
        group=ifelse(original_group_name == ref_group, 'ctrl', 'exp')) %>%
        dplyr::arrange(match(group, c('ctrl', 'exp')), pair)
    if(is.null(ref_group) |
       isFALSE(ref_group %in% unique(group_info$original_group_name))){
        stop("ref_group must be one of the group names in the 'group' column ",
        "of the group information table.")
    }

    ## Paired sample
    paired_sample <- ifelse(all(is.na(group_info$pair)), FALSE, TRUE)
    #### Data transformation ####
    abundance <- .transform(abundance, transform)

    #### Gather the abundance data ####
    sub_abund <- abundance %>% dplyr::filter(feature == .env$feature) %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name')

    #### Two groups comparison ####
    comp.res <- .box_stat_twoGroup(
        sub_abund, group_info, ref_group, paired_sample, test)
    #### Box plot ####
    # if(!is.na(comp.res$res_table$pval)){
        ## Label of y axis
        if(transform == 'log10'){
            Y <- 'Log10-transformed abundance'
        }else if(transform %in% c('cube', 'sqare')){
            Y <- paste0(stringr::str_to_title(transform),
                        '-transformed abundance')
        }else{
            Y <- 'Processed abundance'
        }
        box.res <- .box_plot_twoGroup(
            sub_abund, group_info, res_table=comp.res$res_table,
            test_res=comp.res$test_res, plot_title=unique(sub_abund$feature),
            plot_ylabel=Y)
        return(list(static_boxPlot=box.res$box.plot,
                    table_boxplot=box.res$box.tab,
                    table_stat=comp.res$res_table))
}

.box_stat_twoGroup <- function(sub_abund, group_info, ref_group,
                               paired_sample, test){
    rstatix.ref.group <- setdiff(
        unique(group_info$original_group_name), ref_group)
    if(test == 't-test'){
        test_res <- tryCatch({
            test_res <- rstatix::t_test(
                sub_abund, abund ~ original_group_name, paired=paired_sample,
                ref.group=rstatix.ref.group, var.equal=TRUE)
        }, error=function(e){
            test_res <- NULL
        })
    }else if(test == 'Wilcoxon test'){
        test_res <- tryCatch({
            test_res <- rstatix::wilcox_test(
                sub_abund, exp ~ original_group_name, paired=paired_sample,
                ref.group=rstatix.ref.group)
        }, error=function(e){
            test_res <- NULL
        })
    }
    if(!is.null(test_res)){
        test.table <- test_res %>%
            dplyr::mutate(
                contrast=paste0(group1, '-', group2),
                feature=unique(sub_abund$feature),
                method=test) %>%
            dplyr::select(feature, contrast, statistic, 'pval'='p', method)
    }else{
        test.table <- data.frame(
            feature=unique(sub_abund$feature), contrast=NA,
            statistic=NA, pval=NA, method=test, stringsAsFactors=FALSE)
    }
    return(list(res_table=test.table,
                test_res=test_res))
}

.box_plot_twoGroup <- function(
        sub_abund, group_info, res_table, test_res, plot_title, plot_ylabel){
    if(!is.null(test_res)){
        # Autocompute P-value Positions For Plotting Significance
        pwc <- test_res %>% rstatix::add_xy_position(x='group') %>%
            rstatix::add_significance(p.col='p')
    }
    ## Group as factor
    fac.level <- group_info %>% dplyr::distinct(original_group_name)
    sub_abund$original_group_name <-
        factor(sub_abund$original_group_name,
               levels=fac.level$original_group_name)
    ## Box plot
    BOX <- ggpubr::ggboxplot(
        data=sub_abund, x='original_group_name', y='abund',
        color='original_group_name', palette='Set2',
        title=plot_title, xlab='Group', ylab=plot_ylabel,
        add='jitter', outlier.shape=NA, legend='') +
        ggplot2::theme(axis.title=ggplot2::element_text(size=16),
                       axis.text=ggplot2::element_text(size=14),
                       plot.title=ggplot2::element_text(size=17, hjust=0.5),
                       plot.subtitle=ggplot2::element_text(size=16, hjust=0.5),
                       plot.caption=ggplot2::element_text(size=14)) +
        ggplot2::labs(subtitle=paste0(unique(res_table$method), ', p = ',
                                      format(unique(res_table$pval), digits=3,
                                             scientific=TRUE)))

    if(!is.null(test_res) && !anyNA(res_table$pval)){
        BOX <- BOX +
            ggpubr::stat_pvalue_manual(
                pwc, label='p.signif', hide.ns=TRUE, size=6, bracket.size=0.5) +
            ggplot2::labs(caption=rstatix::get_pwc_label(pwc))
    }
    return(list(box.plot=BOX, box.tab=sub_abund))
}
