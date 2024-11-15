#' @title plot_enrichment_lsea
#' @description This function plots the result of \code{\link{enrichment_lsea}}.
#' @param lsea_res List. A result list of lipid-set enrichment analysis from \code{\link{enrichment_lsea}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @param char_feature Character. A feature selected by users from \code{char}
#' to visualize the specific plot for the selected category of that characteristic.
#' For example, if char is 'class' and char_feature is 'Cer', the resulting
#' plots will display data for 'Cer' within the 'class' category.
#' @return Return a list with 1 static plot and 1 interactive plot.
#' \enumerate{
#' \item static_enrichPlot: a static enriched plot of the selected lipid characteristic and feature.
#' \item interactive_enrichPlot: an interactive enriched plot of the selected lipid characteristic and feature.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se_twoGroup <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' char_list <- list_lipid_char(processed_se_twoGroup)$common_list
#' print(char_list)
#' deSp_se_twoGroup <- deSp_twoGroup(
#'     processed_se_twoGroup, ref_group='ctrl', test='t-test',
#'     significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' lsea_one_twoGroup <- enrichment_lsea(
#'     deSp_se_twoGroup, char='class', rank_by='statistic', significant='pval',
#'     p_cutoff=0.05)
#' lsea_plot_twoGroup <- plot_enrichment_lsea(
#'     lsea_res=lsea_one_twoGroup, char='class', char_feature='TG')
#'
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' char_list <- list_lipid_char(processed_se_multiGroup)$common_list
#' print(char_list)
#' deSp_se_multiGroup <- deSp_multiGroup(
#'     processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
#'     significant='pval', p_cutoff=0.05, transform='log10')
#' lsea_one_multiGroup <- enrichment_lsea(
#'     deSp_se_multiGroup, char='class', rank_by='statistic',
#'     significant='pval', p_cutoff=0.05)
#' lsea_plot_multiGroup <- plot_enrichment_lsea(
#'     lsea_res=lsea_one_multiGroup, char='class', char_feature='Cer')

plot_enrichment_lsea <- function(lsea_res, char, char_feature){
    ## Check parameter
    if(is.null(lsea_res$enrich_result) | is.null(lsea_res$lipid_set)
       | is.null(lsea_res$ranked_list)){
        stop(paste(
            "Incorrect LSEA result object. Please use enrichment_lsea function",
            "to obtain the correct LSEA result object."))
    }
    if(!char %in% lsea_res$enrich_result$characteristic |
       !char_feature %in% lsea_res$enrich_result$charFeature){
        stop(paste(
            "The 'char' and 'char_feature' parameters must be included in the",
            "'enrich_result' returned by enrichment_lsea function."))
    }
    enrich.tab <- lsea_res$enrich_result %>% dplyr::mutate(
        pathway=paste0(characteristic, '_', charFeature), .before=characteristic) %>%
        dplyr::filter(characteristic == char, charFeature == char_feature)
    lipid.set <- lsea_res$lipid_set
    rank.list <- lsea_res$ranked_list
    ## static enrichment plot
    enrich.plot <- fgsea::plotEnrichment(
        lipid.set[[enrich.tab$pathway]], rank.list) +
        ggplot2::labs(title=paste0(
            enrich.tab$characteristic, '\n', enrich.tab$charFeature, '\npval=',
            format(enrich.tab$pval, digits=3, scientific=TRUE), ', padj=',
            format(enrich.tab$padj, digits=3, scientific=TRUE), ', ES=',
            round(enrich.tab$ES, digits=3), ', NES=', round(enrich.tab$NES, digits = 3))) +
        ggplot2::theme(text=ggplot2::element_text(size=15))
    ## interactive enrichment plot
    # enrich.inter <- plotly::ggplotly(enrich.plot) %>%
    #     plotly::layout(title=list(font=list(size=18), y=0.97))
    # data.length <- length(enrich.inter$x$data)
    # for(t in 1:length(enrich.inter$x$data[[data.length]]$text)){
    #     if(!is.na(enrich.inter$x$data[[data.length]]$text[t])){
    #         split.text=rank=lipid=NULL
    #         split.text <- stringr::str_split(
    #             enrich.inter$x$data[[data.length]]$text[t], pattern = '<br />')
    #         rank <- stringr::str_remove(split.text[[1]][1], '.* ') %>%
    #             as.numeric()
    #         lipid <- names(rank.list)[rank]
    #         enrich.inter$x$data[[data.length]]$text[t] <- lipid
    #     }
    # }
    # return(list(static_enrichPlot=enrich.plot,
    #             interactive_enrichPlot=enrich.inter))
    return(static_enrichPlot=enrich.plot)
}
