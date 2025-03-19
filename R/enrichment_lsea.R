#' @title enrichment_lsea
#' @description This function conducts Lipid Set Enrichment Analysis (LSEA), a
#' computational method for determining whether a predefined set of lipids shows
#' statistically significant, concordant differences among two or multiple
#' biological states (e.g., phenotypes).
#' @param deSp_se The resulting SummarizedExperiment object from the differential
#' expression analysis function, such as \code{\link{deSp_twoGroup}} and
#' \code{\link{deSp_multiGroup}}.
#' @param char Character/NULL. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Enter NULL to plot all characteristics.
#' @param rank_by Character. The method to rank lipids. Allowed rankings include
#' "log2FC," "pval," "padj," and "statistic". "log2FC" is only permitted for two-group data.
#' Default is \code{'statistic'}.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. The threshold to distinguish enriched lipid-sets
#' from not-enriched ones. Default is \code{0.05}.
#' @param n_lipid Numeric. The minimum number of lipids in a lipid set to be
#' included in enrichment. Default is \code{2}.
#' @return Return a list of enrichment result, 1 interactive plot, 1 static plot, 1 table, and 2 lists.
#' \enumerate{
#' \item enrich_result: a table of enrichment result.
#' \item static_barPlot: a static bar plot with the top 10 significant up-regulated
#' and down-regulated terms for datasets involving two groups and the top 20 for multiple groups.
#' \item interactive_barPlot: an interactive bar plot with the top 10 significant up-regulated
#' and down-regulated terms for datasets involving two groups and the top 20 for multiple groups.
#' \item table_barPlot: table for plotting bar plots.
#' \item lipid_set: a list of lipid set.
#' \item ranked_list: lipid ranking list.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se_twoGroup <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se_twoGroup <- deSp_twoGroup(
#'     processed_se_twoGroup, ref_group='ctrl', test='t-test',
#'     significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' lsea_all_twoGroup <- enrichment_lsea(
#'     deSp_se_twoGroup, char=NULL, rank_by='statistic', significant='pval',
#'     p_cutoff=0.05)
#' char_list <- list_lipid_char(processed_se_twoGroup)$common_list
#' print(char_list)
#' lsea_one_twoGroup <- enrichment_lsea(
#'     deSp_se_twoGroup, char='class', rank_by='statistic', significant='pval',
#'     p_cutoff=0.05)
#'
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se_multiGroup <- deSp_multiGroup(
#'     processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
#'     significant='pval', p_cutoff=0.05, transform='log10')
#' lsea_all_multiGroup <- enrichment_lsea(
#'     deSp_se_multiGroup, char=NULL, rank_by='statistic', significant='pval',
#'     p_cutoff=0.05)
#' char_list <- list_lipid_char(processed_se_multiGroup)$common_list
#' print(char_list)
#' lsea_one_multiGroup <- enrichment_lsea(
#'     deSp_se_multiGroup, char='class', rank_by='statistic',
#'     significant='pval', p_cutoff=0.05)

enrichment_lsea <- function(deSp_se, char=NULL,
                            rank_by=c('log2FC', 'pval', 'padj', 'statistic'),
                            significant=c('pval', 'padj'), p_cutoff=0.05,
                            n_lipid=2){
    ## Check parameter
    .check_inputSE(deSp_se, metadata_list=c(
        'all_deSp_result', 'sig_deSp_result', 'significant', 'p_cutoff'))
    if(!is.null(char) & isFALSE(char %in% lipidChar$characteristic)){
        stop(paste(
            "The 'char' parameter must be one of the listed lipid",
            "characteristics. Please use the list_lipid_char function to",
            "view which lipid characteristics can be analyzed."))
    }
    if(is.null(rank_by) | isFALSE(rank_by %in% c('log2FC', 'pval', 'padj', 'statistic'))){
        stop("The 'rank_by' parameter must be one of 'log2FC', 'pval', 'padj' or 'statictic'.")
    }
    if(is.null(significant) | isFALSE(significant %in% c('pval', 'padj'))){
        stop("The 'significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1))){
        stop("The 'p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    ## Extract data from SE
    group_info <- .extract_df(deSp_se, type = "group")
    group_type <- .check_nGroup(group_info)
    if(group_type == 'multiple' & rank_by == 'log2FC'){
        stop(paste(
            "Log2 fold change (log2FC) ranking is not suitable for multiple",
            "group data. Please choose p-value (pval), adjusted p-value (padj),",
            "or statistic as the ranking criterion instead."))
    }
    lipid_char <- .extract_df(deSp_se, type = "lipid")
    all_deSp_result <- S4Vectors::metadata(deSp_se)$all_deSp_result

    ## Join lipid_char and all_deSp_result
    deSpecies <- all_deSp_result %>%
        dplyr::left_join(lipid_char, by='feature')

    ## Create gmt file
    gmt <- deSpecies %>%
        dplyr::select(feature, dplyr::all_of(lipidChar$characteristic)) %>%
        tidyr::gather(characteristic, charFeature, -1) %>%
        dplyr::filter(!is.na(charFeature)) %>%
        tidyr::separate_rows(charFeature, sep = '\\|') %>%
        dplyr::mutate(term=paste0(characteristic, '_', charFeature)) %>%
        dplyr::select(term, feature) %>%
        dplyr::group_by(term) %>%
        dplyr::filter(dplyr::n() >= n_lipid) %>%
        dplyr::ungroup()
    if(!is.null(char)) gmt %<>% dplyr::filter(stringr::str_starts(term, paste0(char, '_')))
    lipidset <- split(gmt$feature, gmt$term)
    ## Sort the lipids
    switch(rank_by,
           log2FC=lipid.stat <- deSpecies %>%
               dplyr::select(feature, 'rank'='log2FC') %>%
               dplyr::arrange(dplyr::desc(rank)),
           pval=lipid.stat <- deSpecies %>%
               dplyr::select(feature, 'rank'='negLog10pval') %>%
               dplyr::arrange(dplyr::desc(rank)),
           padj=lipid.stat <- deSpecies %>%
               dplyr::select(feature, 'rank'='negLog10padj') %>%
               dplyr::arrange(dplyr::desc(rank)),
           statistic=lipid.stat <- deSpecies %>%
               dplyr::select(feature, 'rank'='statistic') %>%
               dplyr::arrange(dplyr::desc(rank)))
    lipidrank <- lipid.stat$rank
    names(lipidrank) <- lipid.stat$feature
    ## LSEA
    lsea.res <- suppressWarnings(
        fgsea::fgsea(pathways=lipidset, stats=lipidrank, eps=0, nPermSimple=5000)
    )

    if(significant == 'pval') lsea.res %<>% dplyr::mutate(p=pval)
    if(significant == 'padj') lsea.res %<>% dplyr::mutate(p=padj)
    lsea.res %<>% dplyr::arrange(p) %>%
        tidyr::separate(pathway, into=c('characteristic', 'charFeature'),
                        sep='_', remove=FALSE) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            leadingEdge=paste0(unlist(leadingEdge), collapse = ','),
            significance=ifelse(p < p_cutoff, 'Yes', 'No')) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(lipidChar, by='characteristic')
    if(!is.null(char)){
        ## One characteristic
        bar.res <- .lsea_oneChar(enrich.res=lsea.res, char, rank_by,
                                 significant, p_cutoff)
    }else{
        ## All characteristics
        bar.res <- .lsea_allChar(enrich.res=lsea.res, rank_by,
                                 significant, p_cutoff)
    }
    ## Enrichment table
    enrich.tab <- lsea.res %>% dplyr::arrange(p) %>%
        dplyr::select(significance, aspect, everything(), -pathway, -p, -source)
    return(list(enrich_result=enrich.tab,
                static_barPlot=bar.res$bar.plot,
                interactive_barPlot=bar.res$bar.plotly,
                table_barPlot=bar.res$bar.tab,
                lipid_set=lipidset,
                ranked_list=lipidrank))
}

.lsea_allChar <- function(enrich.res, rank_by, significant, p_cutoff){
    ## Top 10 of all characteristics
    if(rank_by == 'log2FC'){
        top10.up <- enrich.res %>% dplyr::filter(p < p_cutoff, NES > 0) %>%
            dplyr::slice_max(order_by=NES, n=10)
        top10.down <- enrich.res %>% dplyr::filter(p < p_cutoff, NES < 0) %>%
            dplyr::slice_min(order_by=NES, n=10)
        candidate.term <- data.table::rbindlist(list(top10.up, top10.down))
    }else{
        candidate.term <- enrich.res %>%
            dplyr::filter(p < p_cutoff, NES > 0) %>%
            dplyr::slice_max(order_by=NES, n=20)
    }
    # bar plot table
    bar.tab <- candidate.term %>%
        dplyr::mutate(
            yText=ifelse(stringr::str_starts(characteristic, 'Total.')
                         | characteristic %in% c('FA', 'FA.C', 'FA.DB', 'FA.OH'),
                         paste0(characteristic, ': ', charFeature), charFeature)) %>%
        dplyr::mutate(
            yText=stringr::str_wrap(yText, width=30),
            hover=paste0(
                'Aspect: ', aspect, '<br />Characteristic: ', characteristic,
                '<br />Feature: ', charFeature, '<br />-log10(', significant,
                '): ', format(-log10(p), digits=3), '<br />', significant,
                ': ', format(p, digits=3, scientific=TRUE),
                '<br />Significance: ', significance)) %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        dplyr::select(significance, aspect, characteristic, charFeature,
                      p, NES, yText, hover)
    bar.tab$yText <- factor(bar.tab$yText, levels=rev(bar.tab$yText))
    bar <- ggplot2::ggplot(
        bar.tab, ggplot2::aes(x=NES, y=yText, fill=aspect)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept=0, color='#444444') +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(
            breaks=c('Lipid classification', 'Cellular component', 'Function',
                     'Physical or chemical properties', 'Fatty acid properties'),
            values =c("#66C2A5", '#E78AC3', "#FFD92F", "#FC8D62", "#8DA0CB"))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=16),
                       axis.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=12),
                       legend.text=ggplot2::element_text(size=8),
                       legend.position="bottom",
                       plot.title=ggplot2::element_text(size=18)) +
        ggplot2::labs(x='NES', y='')
    bar.plotly <- plotly::ggplotly(bar) %>%
        plotly::layout(legend=list(y=0, x=-0.15, orientation='h'),
                       margin=list(l=70, r=20, b=60, t=40))
    for(i in 1:(length(bar.plotly$x$data)-1)){
        for(j in 1:length(bar.plotly$x$data[[i]]$text)){
            feature <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][2] %>%
                stringr::str_remove('yText: ') %>% stringr::str_wrap(width=30)
            bar.plotly$x$data[[i]]$text[j] <- bar.tab$hover[bar.tab$yText == feature]
        }
    }
    return(list(bar.plot=bar,
                bar.plotly=bar.plotly,
                bar.tab=bar.tab))
}

.lsea_oneChar <- function(enrich.res, char, rank_by, significant, p_cutoff){
    bar.tab <- enrich.res %>%
        tidyr::separate(pathway, into=c('characteristic', 'charFeature'),
                        sep='_', remove=FALSE) %>%
        dplyr::mutate(
            color=ifelse(
                rank_by == 'log2FC' & p < p_cutoff & NES < 0, 'Down',
                ifelse(p < p_cutoff & NES > 0, 'Up', 'NS')),
            significance=ifelse(p < p_cutoff, 'Yes', 'No'),
            yText=stringr::str_wrap(charFeature, width=30),
            hover=paste0(
                'Aspect: ', aspect, '<br />Characteristic: ', characteristic,
                '<br />Feature: ', charFeature, '<br />-log10(', significant,
                '): ', format(-log10(p), digits=3), '<br />', significant,
                ': ', format(p, digits=3, scientific=TRUE),
                '<br />Significance: ', significance)) %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        dplyr::select(significance, aspect, characteristic, charFeature,
                      p, NES, yText, color, hover)
    bar.tab$yText <- factor(bar.tab$yText, levels=rev(bar.tab$yText))
    bar <- ggplot2::ggplot(
        bar.tab, ggplot2::aes(x=NES, y=yText, fill=color)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept=0, color='#444444') +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(
            breaks=c('Up', 'Down', 'NS'),
            values =c("#FF4500", "#4169E1", "#DDDDDD"))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=16),
                       axis.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=12),
                       legend.text=ggplot2::element_text(size=8),
                       legend.position="bottom",
                       plot.title=ggplot2::element_text(size=18)) +
        ggplot2::labs(x='NES', y='', title=char, fill='significance')
    bar.plotly <- plotly::ggplotly(bar) %>%
        plotly::layout(legend=list(y=0, x=-0.15, orientation='h'),
                       margin=list(l=70, r=20, b=60, t=40))
    for(i in 1:(length(bar.plotly$x$data)-1)){
        for(j in 1:length(bar.plotly$x$data[[i]]$text)){
            feature <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][2] %>%
                stringr::str_remove('yText: ') %>% stringr::str_wrap(width=30)
            bar.plotly$x$data[[i]]$text[j] <- bar.tab$hover[which(bar.tab$yText == feature)]
        }
    }
    return(list(bar.plot=bar,
                bar.plotly=bar.plotly,
                bar.tab=bar.tab))
}
