#' @title heatmap_chain_db
#' @description This function generates a heatmap showing the correlation
#' between the double bond and chain length of lipids.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the chain_db list
#' returned by \code{\link{list_lipid_char}}.
#' @param char_feature Character/NULL. A feature selected by users from \code{char}
#' to visualize the specific plot for the selected category of that characteristic.
#' For example, if char is 'class' and char_feature is 'Cer', the resulting
#' plots will display data for 'Cer' within the 'class' category.
#' Set NULL to prevent selecting any feature as char_feature.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param test Character. The method to use for comparing means.
#' Allowed method include "t-test", "Wilcoxon test", "One-way ANOVA", and "Kruskal–Wallis test".
#' "t-test", "Wilcoxon test" are for two-group data, and "One-way ANOVA"
#' and "Kruskal–Wallis test" are for multi-group data. Default is \code{'t-test'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. Significant level. Default is \code{0.05}.
#' @param FC_cutoff Numeric. Significance of the fold-change, which is only
#' applicable for the two-group data. Default is \code{1}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a list of 2 lists.
#' \enumerate{
#' \item total_chain: the result list of total chain.
#' \item each_chain: the result list of fatty acids chain.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se_twoGroup <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' char_list <- list_lipid_char(processed_se_twoGroup)$chain_db_list
#' print(char_list)
#' heatmap_all_twoGroup <- heatmap_chain_db(
#'     processed_se_twoGroup, char='class', char_feature=NULL, ref_group='ctrl',
#'     test='t-test', significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' heatmap_one_twoGroup <- heatmap_chain_db(
#'     processed_se_twoGroup, char='class', char_feature='PC', ref_group='ctrl',
#'     test='t-test', significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#'
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' char_list <- list_lipid_char(processed_se_multiGroup)$chain_db_list
#' print(char_list)
#' heatmap_all_multiGroup <- heatmap_chain_db(
#'     processed_se_multiGroup, char='class', char_feature=NULL, ref_group='ctrl',
#'     test='One-way ANOVA', significant='pval', p_cutoff=0.05, FC_cutoff=NULL,
#'     transform='log10')
#' heatmap_one_multiGroup <- heatmap_chain_db(
#'     processed_se_multiGroup, char='class', char_feature='PC', ref_group='ctrl',
#'     test='One-way ANOVA', significant='pval', p_cutoff=0.05, FC_cutoff=NULL,
#'     transform='log10')

heatmap_chain_db <- function(
        processed_se, char, char_feature=NULL, ref_group,
        test=c('t-test', 'Wilcoxon test', 'One-way ANOVA', 'Kruskal–Wallis test'),
        significant=c('pval', 'padj'), p_cutoff=0.05, FC_cutoff=1,
        transform=c('none', 'log10', 'square', 'cube')){

    #### check parameter ####
    .check_inputSE(processed_se, metadata_list=NULL)
    if(is.null(char) | isFALSE(.check_char(processed_se, char, type='chain_db'))){
        stop(paste("The 'char' parameter is not in the list of lipid",
                   "characteristics. Please execute list_lipid_char() to view",
                   "the valid lipid characteristics you can input."))
    }
    if(is.null(significant) | isFALSE(significant %in% c('pval', 'padj'))){
        stop("The 'significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1))){
        stop("The 'p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    if(is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube'))){
        stop("The transform parameter must be one of 'none', 'log10', 'square' or 'cube'.")
    }

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(processed_se, type="abundance")
    .check_imputation(abundance)
    if(isTRUE(.check_nor_negative(abundance[-1]))) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    ## Group information
    group_info <- .extract_df(processed_se, type="group")
    paired_sample <- ifelse(!all(is.na(group_info$pair)), TRUE, FALSE)
    n_group <- .check_nGroup(group_info)
    if(is.null(ref_group) | isFALSE(ref_group %in% unique(group_info$group))){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    if(n_group == 'two' && (!is.numeric(FC_cutoff) | isFALSE(FC_cutoff >= 1))){
        stop("The minimum value for the 'FC_cutoff' parameter is 1.")
    }
    if(n_group == 'two'
       & any((is.null(test) | !test %in% c('t-test', 'Wilcoxon test'))) ){
        stop(paste("'One-way ANOVA' and 'Kruskal–Wallis test' are for multiple",
                   "group comparisons. Please choose either 't-test' or",
                   "'Wilcoxon test' for two group analysis."))
    }
    if(n_group == 'multiple'
       & (is.null(test) | !test %in% c('One-way ANOVA', 'Kruskal–Wallis test'))){
        stop(paste("'t-test' and 'Wilcoxon test' are for two group comparisons.",
                   "Please choose either 'One-way ANOVA' or 'Kruskal–Wallis test'",
                   "for multiple group analysis."))
    }
    ## Lipid characteristics table
    lipid_char <- .extract_df(processed_se, type="lipid")
    unique_charFeature <- lipid_char %>%
        dplyr::select('charFeature'=dplyr::all_of(char)) %>%
        dplyr::filter(!is.na(charFeature)) %>%
        dplyr::distinct(charFeature) %>% .$charFeature
    if(!is.null(char_feature) && !char_feature %in% unique_charFeature){
        stop("The 'char_feature' parameter must be one of the valid char features.")
    }
    #### Total chain ####
    total_chain <- .chain_db_analysis(
        type='total', abundance, char, char_feature, group_info, lipid_char,
        ref_group, n_group, test, significant, p_cutoff, FC_cutoff,
        paired_sample, transform)
    if(is.null(total_chain)){
        message(paste("This lipid characteristic or characteristic feature does",
                      "not include any lipids for total chain analysis.",
                      "Therefore, the 'total_chain' object is NULL."))
        total_chain <- "This lipid characteristic or characteristic feature does not include any lipids for total chain analysis. Therefore, the 'total_chain' object is NULL."
    }
    #### Each chain ####
    each_chain <- .chain_db_analysis(
        type='each', abundance, char, char_feature, group_info, lipid_char,
        ref_group, n_group, test, significant, p_cutoff, FC_cutoff,
        paired_sample, transform)
    if(is.null(each_chain)){
        message(paste("This lipid characteristic or characteristic feature does",
                      "not include any lipids for fatty acid chain analysis.",
                      "Therefore, the 'each_chain' object is NULL."))
        each_chain <- "This lipid characteristic or characteristic feature does not include any lipids for fatty acid chain analysis. Therefore, the 'each_chain' object is NULL."
    }
    return(list(total_chain=total_chain,
                each_chain=each_chain))

}

.chain_db_analysis <- function(
        type=c('total', 'each'), abundance, char, char_feature=NULL, group_info,
        lipid_char, ref_group, n_group, test, significant, p_cutoff, FC_cutoff,
        paired_sample, transform){
    switch(type,
           total=lipid.char <- lipid_char %>% dplyr::select(
               feature, 'Characteristic'=dplyr::all_of(char), 'FA'='Total.FA'),
           each=lipid.char <- lipid_char %>% dplyr::select(
               feature, 'Characteristic'=dplyr::all_of(char), 'FA'='FA') %>%
               tidyr::separate_rows(FA, sep='\\|'))
    lipid.char %<>%
        tidyr::separate_rows(Characteristic, sep='\\|') %>%
        dplyr::filter(!is.na(Characteristic)) %>%
        dplyr::filter(!is.na(FA))
    if(nrow(lipid.char) == 0) return(NULL)
    ## Chian-DB abundance
    chain_db_abund <- .chain_db_abundance(
        lipid_char=lipid.char, abundance, group_info, char_feature)
    if(nrow(chain_db_abund$char.abund) == 0) return(NULL)
    ## Statistic test
    abund_data <- chain_db_abund$char.abund.mat %>%
        tibble::column_to_rownames(var='FA') %>% as.matrix()
    lipid_char_table <- chain_db_abund$char.abund %>%
        dplyr::filter(FA %in% rownames(abund_data)) %>%
        dplyr::distinct(FA, Chain, DB) %>%
        dplyr::select('feature'='FA', dplyr::everything()) %>%
        dplyr::arrange(match(feature, rownames(abund_data)))
    if(!is.null(char_feature)) lipid_char_table$Characteristic <- char_feature
    group_info <- group_info %>% dplyr::arrange(match(sample_name, colnames(abund_data)))
    chain_db_se <- SummarizedExperiment::SummarizedExperiment(
        assays=abund_data,
        colData=S4Vectors::DataFrame(group_info),
        rowData=S4Vectors::DataFrame(lipid_char_table))
    switch(n_group,
           two=deSp_se <- deSp_twoGroup(
               processed_se=chain_db_se, ref_group, test, significant,
               p_cutoff, FC_cutoff, transform),
           multiple=deSp_se <- deSp_multiGroup(
               processed_se=chain_db_se, ref_group, test, significant,
               p_cutoff, transform))
    ## Plot heatmap
    all_deSp_result <- S4Vectors::metadata(deSp_se)$all_deSp_result
    deSp_lipid_char <- .extract_df(deSp_se, type="lipid")
    all_deSp_result %<>% dplyr::inner_join(deSp_lipid_char, by='feature') %>%
        dplyr::mutate(p.signif=ifelse(pval < p_cutoff, '*', 'ns'),
                      p.adj.signif=ifelse(padj < p_cutoff, '*', 'ns'))
    switch(significant,
           pval=all_deSp_result %<>%
               dplyr::mutate(star=ifelse(p.signif == 'ns', '', p.signif)),
           padj=all_deSp_result %<>%
               dplyr::mutate(star=ifelse(p.adj.signif == 'ns', '', p.adj.signif)))

    switch(n_group,
           two=heatmap_res <- .chain_db_heatmap_twoGroup(
               all_deSp_result, char, char_feature, type, FC_cutoff),
           multiple=heatmap_res <- .chain_db_heatmap_multiGroup(
               all_deSp_result, char, char_feature, type, significant))

    return(list(static_heatmap=heatmap_res$two.char.heatmap,
                table_heatmap=heatmap_res$two.char.tab,
                processed_abundance=S4Vectors::metadata(deSp_se)$processed_abundance,
                transformed_abundance=.extract_df(deSp_se, type="abundance"),
                chain_db_se=chain_db_se))
}

.chain_db_abundance <- function(lipid_char, abundance, group_info, char_feature){
    if(!is.null(char_feature)) lipid_char %<>% dplyr::filter(Characteristic == char_feature)
    char.abund <- lipid_char %>%
        dplyr::filter(!is.na(FA), feature %in% abundance$feature) %>%
        dplyr::mutate(FA=stringr::str_remove(FA, ';.*')) %>%
        tidyr::separate(FA, into=c('Chain', 'DB'), remove=FALSE) %>%
        dplyr::mutate(Chain=as.numeric(Chain),
                      DB=as.numeric(DB)) %>%
        dplyr::left_join(abundance, by='feature') %>%
        tidyr::gather(sample_name, abundance, -1:-5) %>%
        dplyr::left_join(group_info, by='sample_name') %>%
        dplyr::group_by(FA, Chain, DB, sample_name, label_name, group) %>%
        dplyr::summarise(abundance=sum(abundance, na.rm=TRUE), .groups='drop')
    char.abund.mat <- char.abund %>%
        dplyr::select(FA, sample_name, abundance) %>%
        tidyr::spread(sample_name, abundance)
    return(list(char.abund=char.abund,
                char.abund.mat=char.abund.mat))
}

.chain_db_heatmap_twoGroup <- function(
        all_deSp_result, char, char_feature, type, FC_cutoff){

    two.char.tab <- all_deSp_result %>%
        dplyr::mutate(color=ifelse(abs(log2FC) > log2(FC_cutoff), log2FC, NA))

    two.char.heatmap <- two.char.tab %>%
        ggplot2::ggplot(ggplot2::aes(x=Chain, y=DB, width=1, height=1)) +
        ggplot2::geom_tile(ggplot2::aes(fill=color), colour='black') +
        ggplot2::geom_text(ggplot2::aes(label=star), vjust=0.75, size=5) +
        ggplot2::scale_fill_gradient2(
            low="#4169E1", mid="#E5E7E9", high="#FF4500", na.value='#DDDDDD', midpoint=0)+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.text=ggplot2::element_text(size=16),
                       plot.title=ggplot2::element_text(size=22, hjust=0.5),
                       plot.subtitle=ggplot2::element_text(size=20, hjust=0.5),
                       axis.title=ggplot2::element_text(size=18),
                       legend.title=ggplot2::element_text(size=16),
                       legend.text=ggplot2::element_text(size=14))+
        ggplot2::scale_x_continuous(breaks=seq(2, 100, by=2)) +
        ggplot2::scale_y_continuous(breaks=seq(0, 40, by=1)) +
        ggplot2::labs(
            x=ifelse(type == 'each', 'Fatty acid chain length',
                     paste(stringr::str_to_title(type), 'chain length')),
            y=ifelse(type == 'each', 'Fatty acid double bond',
                     paste(stringr::str_to_title(type), 'double bond')),
            fill='Log2FC')
    if(!is.null(char_feature)){
        two.char.heatmap <- two.char.heatmap +
            ggplot2::labs(title=char_feature, subtitle=char)
    }else{
        two.char.heatmap <- two.char.heatmap +
            ggplot2::labs(title=char)
    }
    return(list(two.char.heatmap=two.char.heatmap,
                two.char.tab=two.char.tab))
}

.chain_db_heatmap_multiGroup <- function(
        all_deSp_result, char, char_feature, type, significant){

    switch(significant,
           pval=two.char.tab <- all_deSp_result %>% dplyr::mutate(color=negLog10pval),
           padj=two.char.tab <- all_deSp_result %>% dplyr::mutate(color=negLog10padj))

    two.char.heatmap <- two.char.tab %>%
        ggplot2::ggplot(ggplot2::aes(x=Chain, y=DB, width=1, height=1)) +
        ggplot2::geom_tile(ggplot2::aes(fill=color), colour='black') +
        ggplot2::geom_text(ggplot2::aes(label=star), vjust=0.75, size=5) +
        ggplot2::scale_fill_gradient2(high="#FF4500", na.value='gray')+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.text=ggplot2::element_text(size=16),
                       plot.title=ggplot2::element_text(size=22, hjust=0.5),
                       plot.subtitle=ggplot2::element_text(size=20, hjust=0.5),
                       axis.title=ggplot2::element_text(size=18),
                       legend.title=ggplot2::element_text(size=16),
                       legend.text=ggplot2::element_text(size=14))+
        ggplot2::scale_x_continuous(breaks=seq(2, 100, by=2)) +
        ggplot2::scale_y_continuous(breaks=seq(0, 40, by=1)) +
        ggplot2::labs(
            x=ifelse(type == 'each', 'Fatty acid chain length',
                     paste(stringr::str_to_title(type), 'chain length')),
            y=ifelse(type == 'each', 'Fatty acid double bond',
                     paste(stringr::str_to_title(type), 'double bond')),
            fill=paste0('-log10(', significant, ')'))
    if(!is.null(char_feature)){
        two.char.heatmap <- two.char.heatmap +
            ggplot2::labs(title=char_feature, subtitle=char)
    }else{
        two.char.heatmap <- two.char.heatmap +
            ggplot2::labs(title=char)
    }
    return(list(two.char.heatmap=two.char.heatmap,
                two.char.tab=two.char.tab))
}
