#' @title enrichment_ora
#' @description This function conducts Over-Representation Analysis (ORA) to
#' determine whether significant lipid species are enriched in specific lipid class categories.
#' @param deSp_se The resulting SummarizedExperiment object from the differential
#' expression analysis function, such as \code{\link{deSp_twoGroup}} and
#' \code{\link{deSp_multiGroup}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}. Enter NULL to plot all characteristics.
#' @param significant Character. The p-value to be used for the statistically
#' significant. Must be one of "pval" or "padj". Default is \code{'pval'}.
#' @param p_cutoff Numeric. The threshold to distinguish enriched lipid-sets
#' from not-enriched ones. Default is \code{0.05}.
#' @return Return a list of enrichment result, 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item enrich_result: a table of enrichment result.
#' \item static_barPlot: a static bar plot that classifies significant lipid species into 'up-regulated' or 'down-regulated' categories.
#' \item interactive_barPlot: an interactive bar plot that classifies significant lipid species into 'up-regulated' or 'down-regulated' categories.
#' \item table_barPlot: table for plotting bar plots.
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
#' ora_all_twoGroup <- enrichment_ora(
#'     deSp_se_twoGroup, char=NULL, significant='pval', p_cutoff=0.05)
#' char_list <- list_lipid_char(processed_se_twoGroup)$common_list
#' print(char_list)
#' ora_one_twoGroup <- enrichment_ora(
#'     deSp_se_twoGroup, char='class', significant='pval', p_cutoff=0.05)
#'
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se_multiGroup <- deSp_multiGroup(
#'     processed_se_multiGroup, ref_group='ctrl', test='One-way ANOVA',
#'     significant='pval', p_cutoff=0.05, transform='log10')
#' ora_all_multiGroup <- enrichment_ora(
#'     deSp_se_multiGroup, char=NULL, significant='pval', p_cutoff=0.05)
#' char_list <- list_lipid_char(processed_se_multiGroup)$common_list
#' print(char_list)
#' ora_one_multiGroup <- enrichment_ora(
#'     deSp_se_multiGroup, char='class', significant='pval', p_cutoff=0.05)

enrichment_ora <- function(deSp_se, char=NULL, significant=c('pval', 'padj'),
                           p_cutoff=0.05){
    ## Check parameter
    .check_de_outputSE(deSp_se, de_type="deSp")
    if(!is.null(char) & isFALSE(char %in% lipidChar$characteristic)){
        stop("The 'char' parameter must be one of the listed lipid ",
             "characteristics. Please use the list_lipid_char function to ",
             "view which lipid characteristics can be analyzed.")
    }
    if(is.null(significant) | isFALSE(significant %in% c('pval', 'padj'))){
        stop("The 'significant' parameter must be either 'pval' or 'padj'.")
    }
    if(!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1))){
        stop("The 'p_cutoff' parameter must be a numeric value between 0 and 1.")
    }
    ## Extract data from SE
    lipid_char <- .extract_df(deSp_se, type = "lipid")
    all_deSp_result <- S4Vectors::metadata(deSp_se)$all_deSp_result
    sig_deSp_result <- S4Vectors::metadata(deSp_se)$sig_deSp_result
    de_significant <- S4Vectors::metadata(deSp_se)$significant
    de_p_cutoff <- S4Vectors::metadata(deSp_se)$p_cutoff
    de_FC_cutoff <- S4Vectors::metadata(deSp_se)$FC_cutoff
    ## Join lipid_char and all_deSp_result
    deSpecies <- all_deSp_result %>%
        dplyr::left_join(lipid_char, by='feature') %>%
        dplyr::mutate(direction=NA, .before='method')
    if(de_significant == 'pval') deSpecies %<>% dplyr::mutate(p=pval)
    if(de_significant == 'padj') deSpecies %<>% dplyr::mutate(p=padj)
    ## Subset significant lipids
    if(!is.null(de_FC_cutoff)){
        deSpecies %<>% dplyr::mutate(direction=ifelse(log2FC > 0, 'Up', 'Down'))
        sigLipid <- deSpecies %>% dplyr::filter(
            p < de_p_cutoff, abs(log2FC) > log2(de_FC_cutoff))
    }else{
        sigLipid <- deSpecies %>% dplyr::filter(p < de_p_cutoff)
    }
    ## Remove NA and separate rows
    all.lipid <- deSpecies %>%
        dplyr::select(feature, direction, dplyr::all_of(lipidChar$characteristic)) %>%
        tidyr::gather(characteristic, charFeature, -seq_len(2)) %>%
        dplyr::filter(!is.na(charFeature)) %>%
        tidyr::separate_rows(charFeature, sep = '\\|')
    sig.lipid <- all.lipid %>% dplyr::filter(feature %in% sigLipid$feature)
    if(nrow(sig.lipid) == 0){
        return(warning(
            'Your data does not contain any significant lipids. Please reset ',
            'the cutoff or review the differential expression results.'))
    }
    ## Fisher's exact test
    fisher.res <- .lipoFisher(sig.lipid, all.lipid)
    if(significant == 'pval') fisher.res %<>% dplyr::mutate(p=pval)
    if(significant == 'padj') fisher.res %<>% dplyr::mutate(p=padj)
    fisher.res %<>%
        dplyr::mutate(significance=ifelse(p < p_cutoff, 'Yes', 'No')) %>%
        dplyr::left_join(lipidChar, by='characteristic')
    if(!is.null(char)){
        ## One characteristic
        fisher.res %<>% dplyr::filter(characteristic == char)
        bar.res <- .ora_oneChar(enrich.res=fisher.res, char, significant)
    }else{
        ## All characteristics
        bar.res <- .ora_allChar(enrich.res=fisher.res, significant)
    }
    ## Enrichment table
    enrich.tab <- fisher.res %>% dplyr::arrange(p) %>%
        dplyr::select(significance, aspect, characteristic, charFeature,
                      'condition'='direction', pval, padj, negLog10pval,
                      negLog10padj, InCharDeLipid, notInCharDeLipid,
                      InCharLipid, notInCharLipid, Lipids)
    if(is.null(de_FC_cutoff)) enrich.tab %<>% dplyr::select(-condition)
    return(list(enrich_result=enrich.tab,
                static_barPlot=bar.res$bar.plot,
                interactive_barPlot=bar.res$bar.plotly,
                table_barPlot=bar.res$bar.tab))
}

.lipoFisher <- function(sig.lipid, all.lipid){
    # All lipids in the one class
    all.class <- all.lipid %>%
        dplyr::group_by(characteristic, charFeature, direction)%>%
        dplyr::summarise(nLipid=dplyr::n(), .groups='drop')
    # Significant lipids in the one class
    sig.class <- sig.lipid %>%
        dplyr::arrange(charFeature) %>%
        dplyr::group_by(characteristic, charFeature, direction)%>%
        dplyr::summarise(InCharDeLipid=dplyr::n(),
                         Lipids=paste0(feature, collapse=','), .groups='drop')%>%
        dplyr::left_join(all.class, by=c('characteristic', 'charFeature', 'direction'))
    # Number of all significant lipids & all lipids
    n.lipid <- sig.class %>%
        dplyr::group_by(characteristic) %>%
        dplyr::summarise(all.sig = sum(InCharDeLipid, na.rm=TRUE),
                         all.lipid = sum(nLipid, na.rm=TRUE), .groups='drop')
    # Calculate the 2*2 contingency table
    sig.class %<>% dplyr::left_join(n.lipid, by='characteristic') %>%
        dplyr::mutate(notInCharDeLipid=all.sig-InCharDeLipid,
                      InCharLipid=nLipid-InCharDeLipid,
                      notInCharLipid=all.lipid-InCharDeLipid-notInCharDeLipid-InCharLipid,
                      pval=NA)
    # Fisher's exact test
    for(i in seq_len(nrow(sig.class))){
        sig.class$pval[i] <- stats::fisher.test(matrix(c(sig.class$InCharDeLipid[i],
                                                      sig.class$notInCharDeLipid[i],
                                                      sig.class$InCharLipid[i],
                                                      sig.class$notInCharLipid[i]),
                                                    nrow=2),
                                             alternative = 'greater')$p.value
    }
    sig.class %<>% dplyr::group_by(characteristic) %>%
        dplyr::mutate(padj=p.adjust(pval, method='BH')) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(negLog10pval=(-log10(pval)),
                      negLog10padj=(-log10(padj)))
    return(sig.class)
}

.ora_allChar <- function(enrich.res, significant){
    if(all(is.na(enrich.res$direction))){
        plot.tab <- enrich.res %>% dplyr::filter(significance == 'Yes') %>%
            dplyr::arrange(p) %>% dplyr::top_n(n=-20, wt=p)
    }else{
        plot.tab <- enrich.res %>% dplyr::filter(significance == 'Yes') %>%
            dplyr::arrange(p) %>% dplyr::group_by(direction) %>%
            dplyr::top_n(n=-10, wt=p)
    }
    if(nrow(plot.tab) == 0) return(NULL)
    x.max <- max(ceiling(plot.tab$negLog10pval))
    if(x.max == 0){
        x.label <- 0
    }else{
        x.label <- as.character(c(rev(seq_len(x.max)), 0, seq_len(x.max)))
    }
    plot.tab %<>% dplyr::mutate(
        negLog10p=ifelse(!is.na(direction) & direction == 'Down', log10(p), -log10(p)),
        hover=paste0(
            'Aspect: ', aspect, '<br />Characteristic: ', characteristic,
            '<br />Feature: ', charFeature, '<br />-log10(', significant, '): ',
            format(-log10(p), digits=3), '<br />', significant, ': ',
            format(p, digits=3, scientific=TRUE), '<br />Significance: ', significance),
        yText=ifelse(stringr::str_starts(characteristic, 'Total.')
                     | characteristic %in% c('FA', 'FA.C', 'FA.DB', 'FA.OH'),
                     paste0(characteristic, ': ', charFeature), charFeature)) %>%
        dplyr::mutate(yText=stringr::str_wrap(yText, width=30)) %>%
        dplyr::group_by(characteristic, charFeature) %>%
        dplyr::mutate(yOrder=negLog10p[which(abs(negLog10p) == max(abs(negLog10p)))]) %>%
        dplyr::select(significance, aspect, characteristic, charFeature,
                      direction, p, negLog10p, yText, yOrder, hover)
    bar <- ggplot2::ggplot(
        plot.tab, ggplot2::aes(x=negLog10p, y=reorder(yText, yOrder, max), fill=aspect)) +
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
        ggplot2::labs(x=paste0('-log10(', significant, ')'), y='') +
        ggplot2::scale_x_continuous(breaks=-x.max:x.max, labels=x.label)
    bar.plotly <- plotly::ggplotly(bar) %>%
        plotly::layout(legend=list(y=0, x=-0.15, orientation='h'),
                       margin=list(l=70, r=20, b=60, t=40))
    for(i in seq_len(length(bar.plotly$x$data)-1)){
        for(j in seq_len(length(bar.plotly$x$data[[i]]$text))){
            logP <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][1] %>%
                stringr::str_remove_all('.* ') %>% as.numeric() %>% round(digits=3)
            feature <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][2] %>%
                stringr::str_remove('.*\\): ') %>% stringr::str_wrap(width=30)
            if(sum(round(plot.tab$negLog10p, digits=3) == logP & plot.tab$yText == feature) == 1 ){
                bar.plotly$x$data[[i]]$text[j] <- plot.tab$hover[which(round(plot.tab$negLog10p, digits=3) == logP
                    & plot.tab$yText == feature)]
            }
        }
    }
    return(list(bar.plot=bar,
                bar.plotly=bar.plotly,
                bar.tab=plot.tab))
}

.ora_oneChar <- function(enrich.res, char, significant){
    if(all(is.na(enrich.res$direction))){
        plot.tab <- enrich.res %>% dplyr::mutate(color=significance)
    }else{
        plot.tab <- enrich.res %>%
            dplyr::mutate(
                color=ifelse(
                    significance == 'Yes' & direction == 'Up', 'Up',
                    ifelse(significance == 'Yes' & direction == 'Down', 'Down', 'NS')))
    }
    x.max <- max(ceiling(plot.tab$negLog10pval))
    if(x.max == 0){
        x.label <- 0
    }else{
        x.label <- as.character(c(rev(seq_len(x.max)), 0, seq_len(x.max)))
    }
    plot.tab %<>% dplyr::mutate(
        negLog10p=ifelse(!is.na(direction) & direction == 'Down', log10(p), -log10(p)),
        hover=paste0(
            'Aspect: ', aspect, '<br />Characteristic: ', characteristic,
            '<br />Feature: ', charFeature, '<br />-log10(', significant, '): ',
            format(-log10(p), digits=3), '<br />', significant, ': ',
            format(p, digits=3, scientific=TRUE), '<br />Significance: ', significance),
        yText=stringr::str_wrap(charFeature, width=30)) %>%
        dplyr::group_by(charFeature) %>%
        dplyr::mutate(yOrder=negLog10p[which(abs(negLog10p) == max(abs(negLog10p)))]) %>%
        dplyr::select(significance, aspect, characteristic, charFeature,
                      direction, p, negLog10p, yText, yOrder, color, hover)
    bar <- ggplot2::ggplot(
        plot.tab, ggplot2::aes(x=negLog10p, y=reorder(yText, yOrder, max), fill=color)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept=0, color='#444444') +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(
            breaks=c('Up', 'Down', 'NS', 'Yes', 'No'),
            values =c("#FF4500", "#4169E1", "#DDDDDD", "#FF4500", "#DDDDDD"))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=16),
                       axis.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=12),
                       legend.text=ggplot2::element_text(size=8),
                       legend.position="bottom",
                       plot.title=ggplot2::element_text(size=18)) +
        ggplot2::labs(x=paste0('-log10(', significant, ')'), y='',
                      title=char, fill='significance') +
        ggplot2::scale_x_continuous(breaks=-x.max:x.max, labels=x.label)
    bar.plotly <- plotly::ggplotly(bar) %>%
        plotly::layout(legend=list(y=0, x=-0.15, orientation='h'),
                       margin=list(l=70, r=20, b=60, t=40))
    for(i in seq_len(length(bar.plotly$x$data)-1)){
        for(j in seq_len(length(bar.plotly$x$data[[i]]$text))){
            # logP <- stringr::str_split(
            #     bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][1] %>%
            #     stringr::str_remove_all('.* ') %>% as.numeric() %>% round(digits=3)
            feature <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][2] %>%
                stringr::str_remove('.*\\): ') %>% stringr::str_wrap(width=30)
            colors <- stringr::str_split(
                bar.plotly$x$data[[i]]$text[j], pattern='<br />', n=3)[[1]][3] %>%
                stringr::str_remove('.*: ')
            if(sum(plot.tab$color == colors & plot.tab$yText == feature) == 1 ){
                bar.plotly$x$data[[i]]$text[j] <- plot.tab$hover[which(plot.tab$color == colors
                    & plot.tab$yText == feature)]
            }
        }
    }
    return(list(bar.plot=bar,
                bar.plotly=bar.plotly,
                bar.tab=plot.tab))
}



