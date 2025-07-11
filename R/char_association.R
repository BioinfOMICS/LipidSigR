#' @title char_association
#' @description This function shows the significant lipid species based on
#' different lipid characteristics and visualizes the difference between control
#' and experimental groups.
#' @param deSp_se The resulting SummarizedExperiment object from the
#' differential expression analysis function, such as
#' \code{\link{deSp_twoGroup}} and \code{\link{deSp_multiGroup}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @return Return a list with 3 tables, 3 interactive plots, and 3 static plots.
#' \enumerate{
#' \item interactive_barPlot & static_barPlot: a bar plot distinguishes
#' significant groups (values) exhibiting a mean fold change greater than 1.
#' (NOTE: The bar chart is not generate for multiple-group data.)
#' \item interactive_lollipop & static_lollipop: a lollipop plot compares all
#' significant groups within the selected characteristic by log2(fold change)
#' for two-group data and -log10(p-value) for multiple-group data.
#' \item interactive_wordCloud & static_wordCloud: a word cloud visualizes the
#' frequency of each group (value) associated with the chosen characteristic.
#' \item table_barPlot: table for plotting bar plots.
#' \item table_lollipop: table for plotting lollipop plots.
#' \item table_wordCloud: table for plotting word cloud.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10' )
#' deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
#'     significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' result <- char_association(deSp_se, char='class')
char_association <- function(deSp_se, char){
    ## check data
    .check_de_outputSE(deSp_se, de_type="deSp")
    if (is.null(char) | isFALSE(.check_char(deSp_se, char, type='common'))) {
        stop("Wrong char input, you can view the available char list by list_lipid_char function.")
    }

    # Subset the data
    abundance_data <- .extract_df(deSp_se, type="abundance")
    lipidChar <- .extract_df(deSp_se, type="lipid")
    group_info <- .extract_df(deSp_se, type="group")
    significant <- S4Vectors::metadata(deSp_se)[["significant"]]
    FC_cutoff <- S4Vectors::metadata(deSp_se)[["FC_cutoff"]]
    sig_table <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]]
    sigDeSpecies <- sig_table%>% dplyr::left_join(lipidChar, by='feature')

    n_group <- ifelse(length(unique(group_info$group))==2, 'two', 'multiple')
    switch(n_group,
           two=plot.tab <- sigDeSpecies %>%
               dplyr::select(
                   feature, FC, log2FC, 'characteristic'=all_of(char),
                   pval, negLog10pval, padj, negLog10padj),
           multiple=plot.tab <- sigDeSpecies %>%
               dplyr::select(
                   feature, 'characteristic'=all_of(char), pval, negLog10pval,
                   padj, negLog10padj))
    plot.tab %<>% dplyr::filter(!is.na(characteristic)) %>%
        tidyr::separate_rows(characteristic, sep='\\|')
    if (nrow(plot.tab) == 0) {
        stop("Insufficient number of significant lipids.")
    }
    plot.tab$characteristic <- as.character(plot.tab$characteristic)
    # Plot
    if(n_group == 'two'){
        # Bar plot
        bar.res <- .sigDeSpeciesBar(plot.tab, char, FC_cutoff)
        interactive_barPlot <- bar.res$bar.plotly
        static_barPlot <- bar.res$bar.plot
        table_barPlot <- bar.res$bar.tab
        # Lollipop plot
        lolli.res <- .sigDeSpeciesLollipopTwoGroups(plot.tab, char)
    }else{
        interactive_barPlot <-
            "Only provided when two-group differential analysis data is input."
        static_barPlot <-
            "Only provided when two-group differential analysis data is input."
        table_barPlot <-
            "Only provided when two-group differential analysis data is input."
        # Lollipop plot
        lolli.res <- .sigDeSpeciesLollipopMultiGroups(plot.tab, char, significant)
    }
    # Word cloud
    wc.res <- .sigDeSpeciesWordCloud(plot.tab)
    return(list(
        interactive_barPlot=interactive_barPlot,
        interactive_lollipop=lolli.res$lolli.plotly,
        interactive_wordCloud=wc.res$wc,
        static_barPlot=static_barPlot, static_lollipop=lolli.res$lolli.plot,
        static_wordCloud=wc.res$wcStatic,
        table_barPlot=table_barPlot, table_lollipop=lolli.res$lolli.tab,
        table_wordCloud=wc.res$wc.tab))
} #Function: Sig_lipid_feature

.sigDeSpeciesBar <- function(plot.tab, char, FC_cutoff){
    # Table for bar plot
    bar.tab <- plot.tab %>% dplyr::group_by(characteristic) %>%
        dplyr::summarise(log2FC.mean=mean(log2FC, na.rm=TRUE),
                  log2FC.sd=sd(log2FC, na.rm=TRUE),
                  log2FC.direction=log2FC.mean/abs(log2FC.mean)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(significant=ifelse(abs(log2FC.mean) > log2(FC_cutoff), 'Yes', 'No'),
               log2FC.sd=ifelse(is.na(log2FC.sd), 0, log2FC.sd),
               hover=paste0(
                   'Characteristic: ', char, '<br />Feature: ', characteristic,
                   '<br />Mean of log2FC: ',
                   format(log2FC.mean, digits=3, scientific=TRUE),
                   '<br />Significance: ', significant),
               x=stringr::str_wrap(characteristic, width=10)) %>%
        dplyr::arrange(desc(log2FC.mean))
    bar.tab$x <- factor(bar.tab$x, levels=bar.tab$x)
    # Bar plot
    bar <- ggplot2::ggplot(bar.tab,
                           ggplot2::aes(x=x, y=log2FC.mean, fill=significant)) +
        ggplot2::geom_col() +
        ggplot2::geom_errorbar(
            ggplot2::aes(min=log2FC.mean,
                         max=log2FC.mean+(log2FC.sd*log2FC.direction))) +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(breaks=c('Yes', 'No'),
                                   values=c('#FF4500', '#666666')) +
        ggplot2::labs(x='', y='log2(Fold Change)',
                      fill=paste0('> ', FC_cutoff, ' FC')) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=16, hjust=0.5),
                       axis.title=ggplot2::element_text(size=14),
                       axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1, size=12),
                       axis.text.y=ggplot2::element_text(size=12))
    # Plotly
    bar.plotly <- plotly::ggplotly(bar)
    # for(t in 1:length(bar.plotly$x$data)){
    #     for(l in 1:length(bar.plotly$x$data[[t]]$text)){
    #         c <- stringr::str_split(bar.plotly$x$data[[t]]$text[l], '<br />', 2)[[1]][1] %>% stringr::str_remove('x: ')
    #         f <- stringr::str_split(bar.plotly$x$data[[t]]$text[l], '<br />', 2)[[1]][1] %>% stringr::str_remove('x: ')
    #         bar.plotly$x$data[[t]]$text[l] <- bar.tab$hover[which(bar.tab$x == c)]
    #     }
    # }
    for(t in seq_len(length(bar.plotly$x$data))){
        for(l in seq_len(length(bar.plotly$x$data[[t]]$text))){
            bar.plotly$x$data[[t]]$text[l] <- bar.plotly$x$data[[t]]$text[l] %>%
                # First remove the log2FC.mean line
                gsub("log2FC\\.mean\\s*\\+\\s*\\(log2FC\\.sd\\s*\\*\\s*log2FC\\.direction\\):\\s*-?[0-9.]+<br />", "", .)
            if (isFALSE(startsWith(bar.plotly$x$data[[t]]$text[l], "x: ")) ) {
                bar.plotly$x$data[[t]]$text[l] <- stringr::str_split(bar.plotly$x$data[[t]]$text[l], '<br />', 2)[[1]][-1]
            }
            bar.plotly$x$data[[t]]$text[l] <- gsub('x: ', 'Characteristic: ', bar.plotly$x$data[[t]]$text[l])
            bar.plotly$x$data[[t]]$text[l] <- gsub('log2FC.mean: ','mean(log2FC) :', bar.plotly$x$data[[t]]$text[l])
        }
    }
    return(list(bar.plotly=bar.plotly, bar.plot=bar, bar.tab=bar.tab))
}

.sigDeSpeciesLollipopTwoGroups <- function(plot.tab, char){
    # Table for lollipop plot
    lolli.tab <- plot.tab %>%
        dplyr::mutate(hover=paste0(
            "Characteristic: ", char, '<br />Feature: ', feature,
            '<br />log2FC : ', format(log2FC, digits=3)),
            x=stringr::str_wrap(characteristic, width=10)) %>%
        dplyr::arrange(dplyr::desc(log2FC))
    lolli.tab$x <- factor(lolli.tab$x, levels=rev(unique(lolli.tab$x)))
    # Lollipop plot
    lollipop <- ggpubr::ggdotchart(
        lolli.tab, combine=TRUE, x="x", y="log2FC", rotate=TRUE, color='characteristic',
        sorting="none", add="segments", dot.size=5, legend.title=char, xlab='',
        ylab="log2(Fold change)", legend="right", ggtheme=ggpubr::theme_pubr())
    # Lollipop plotly
    lollipop.plotly <- plotly::ggplotly(lollipop) %>% plotly::layout(showlegend=FALSE)
    for(t in seq_len(length(lollipop.plotly$x$data))){
        if(t == 1){
            lollipop.plotly$x$data[[t]]$text <- NULL
        }else{
            for(l in seq_len(length(lollipop.plotly$x$data[[t]]$text))){
                c <- stringr::str_split(lollipop.plotly$x$data[[t]]$text[l], '<br />', 2)[[1]][1] %>%
                    stringr::str_remove('.*: ')
                lollipop.plotly$x$data[[t]]$text[l] <- lolli.tab$hover[which(lolli.tab$x == c)][l]
            }
        }
    }
    return(list(lolli.plotly=lollipop.plotly, lolli.plot=lollipop, lolli.tab=lolli.tab))
}

.sigDeSpeciesLollipopMultiGroups <- function(plot.tab, char, significant){
    # Table for lollipop plot
    if(significant == 'pval'){
        lolli.tab <- plot.tab %>%
            dplyr::mutate(
                y=negLog10pval,
                hover=paste0("Characteristic: ", char, '<br />Feature: ', feature,
                                '<br />p-value : ', format(pval, digits=3, scientific=TRUE),
                                '<br />-log10pval : ', format(negLog10pval, digits=3)),
                          x=stringr::str_wrap(characteristic, width=10)) %>%
            dplyr::arrange(dplyr::desc(negLog10pval))
    }else if(significant == 'padj'){
        lolli.tab <- plot.tab %>%
            dplyr::mutate(
                y=negLog10padj,
                hover=paste0("Characteristic: ", char, '<br />Feature: ', characteristic,
                               '<br />Adjusted p-value : ', format(padj, digits=3, scientific=TRUE),
                               '<br />-log10padj : ', format(negLog10padj, digits=3)),
                          x=stringr::str_wrap(characteristic, width=10)) %>%
            dplyr::arrange(dplyr::desc(negLog10padj))
    }
    lolli.tab$x <- factor(lolli.tab$x, levels=rev(unique(lolli.tab$x)))
    # Lollipop plot
    lollipop <- ggpubr::ggdotchart(
        lolli.tab, combine=TRUE, x="x", y='y', rotate=TRUE, color='characteristic',
        sorting="none", add="segments", dot.size=5, legend.title=char, xlab='',
        ylab=paste0("-log10(", significant, ")"), legend="right", ggtheme=ggpubr::theme_pubr())
    # Lollipop plotly
    lollipop.plotly <- plotly::ggplotly(lollipop)
    for(t in seq_len(length(lollipop.plotly$x$data))){
        if(t == 1){
            lollipop.plotly$x$data[[t]]$text <- NULL
        }else{
            for(l in seq_len(length(lollipop.plotly$x$data[[t]]$text))){
                c <- stringr::str_split(lollipop.plotly$x$data[[t]]$text[l], '<br />', 2)[[1]][1] %>%
                    stringr::str_remove('.*: ')
                lollipop.plotly$x$data[[t]]$text[l] <- lolli.tab$hover[which(lolli.tab$x == c)][l]
            }
        }
    }
    return(list(lolli.plotly=lollipop.plotly, lolli.plot=lollipop, lolli.tab=lolli.tab))
}

.sigDeSpeciesWordCloud <- function(plot.tab){
    wc.tab <- plot.tab %>% dplyr::group_by(characteristic) %>%
        dplyr::summarise(freqs=dplyr::n())
    if(nrow(wc.tab)==1){
        wc <- hwordcloud::hwordcloud(text=c(wc.tab$characteristic, NA),
                                     size=c(wc.tab$freqs, NA), theme="sunset")
    }else{
        wc <- hwordcloud::hwordcloud(text=wc.tab$characteristic, size=wc.tab$freqs, theme="sunset")
    }

    wordcloud::wordcloud(
        wc.tab$characteristic, wc.tab$freqs, min.freq=1, random.order=FALSE,
        ordered.colors=FALSE, colors=grDevices::rainbow(nrow(wc.tab)))
    wcStatic <- grDevices::recordPlot()
    grDevices::dev.off()
    return(list(wc=wc, wcStatic=wcStatic, wc.tab=wc.tab))
}
