#' @title plot_deSp_multiGroup
#' @description This function is for plotting the results of lipid species differential expression analysis.
#' @param deSp_se A SummarizedExperiment object with results computed by \code{\link{deSp_multiGroup}}.
#' @return Return a list of 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item interactive_de_lipid & static_de_lipid: a lollipop chart reveals the top 20 lipid species that pass chosen cut-off.
#' \item interactive_dotPlot & static_dotPlot: a dot plot reveals the all lipid species that pass chosen cut-off.
#' \item table_de_lipid: table for plotting DE lipid plot.
#' \item table_dotPlot: table for plotting dot plot.
#' }
#' @export
#' @examples
#' data("se_multiGroup")
#' processed_se <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se <- deSp_multiGroup(processed_se, ref_group='ctrl', test='One-way ANOVA',
#'     significant='pval', p_cutoff=0.05, transform='log10')
#' deSp_plot <- plot_deSp_multiGroup(deSp_se)

plot_deSp_multiGroup <- function(deSp_se){
    .check_de_outputSE(deSp_se, de_type="deSp")
    ## Extract data from SE
    lipid_char <- .extract_df(deSp_se, type = "lipid")
    all_deSp_result <- S4Vectors::metadata(deSp_se)$all_deSp_result
    sig_deSp_result <- S4Vectors::metadata(deSp_se)$sig_deSp_result
    significant <- S4Vectors::metadata(deSp_se)$significant
    ## Check results
    if(is.character(sig_deSp_result)){
        warning("This case does not include any significant lipids; ",
                "therefore, all lipids displayed on the plot are ",
                "non-significant.")
    }
    if(nrow(all_deSp_result) > 0){
        ## Lollipop plot
        lollipop <- .lollipop_multiGroup(all_deSp_result, significant)
        ## Dot plot
        dot <- .dot_multiGroup(all_deSp_result, lipid_char, significant)
        return(list(static_de_lipid=lollipop$static_de_lipid,
                    static_dotPlot=dot$static_dotPlot,
                    interactive_de_lipid=lollipop$interactive_de_lipid,
                    interactive_dotPlot=dot$interactive_dotPlot,
                    table_de_lipid=lollipop$table_de_lipid,
                    table_dotPlot=dot$table_dotPlot))
    }else{
        message("This case does not include any lipids for generating ",
                "plots; therefore, the object is NULL.")
        return(NULL)
    }
}


.lollipop_multiGroup <- function(deSp_result, significant=c('pval', 'padj')){
    if(significant == 'pval'){
        top.tab <- deSp_result %>% as.data.frame() %>%
            dplyr::top_n(n=20, wt=negLog10pval) %>%
            dplyr::mutate(
                hover=paste0('Lipid: ', feature, '<br />-log10(pval): ',
                             format(negLog10pval, digits=3), '<br />pval: ',
                             format(pval, digits=3, scientific=TRUE)))
        lolli.plot <- ggpubr::ggdotchart(
            top.tab, x="feature", y="negLog10pval", rotate=TRUE, color="negLog10pval",
            sorting="descending", add="segments", dot.size=4, legend="right",
            legend.title="-log10(pval)", xlab=" ", ylab="-log10(pval)",
            ggtheme=ggpubr::theme_pubr())
    }else if(significant == 'padj'){
        top.tab <- deSp_result %>% as.data.frame() %>%
            dplyr::top_n(n=20, wt=negLog10padj) %>%
            dplyr::mutate(
                hover=paste0('Lipid: ', feature, '<br />-log10(padj): ',
                             format(negLog10padj, digits=3), '<br />padj: ',
                             format(padj, digits=3, scientific=TRUE)))
        lolli.plot <- ggpubr::ggdotchart(
            top.tab, x="feature", y="negLog10padj", rotate=TRUE,
            color="negLog10padj", sorting="descending", add="segments",
            dot.size=4, legend="right", legend.title="-log10(padj)", xlab=" ",
            ylab="-log10(padj)", ggtheme=ggpubr::theme_pubr())
    }
    lolli.plot <- lolli.plot +
        ggplot2::scale_colour_gradient(
            low="#c6dbef", high="#08519c", na.value=NA) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=14),
                       axis.title=ggplot2::element_text(size=16),
                       legend.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=14))
    int.lolli <- plotly::ggplotly(lolli.plot)
    n <- t <- NULL
    for(n in seq_len(2)){
        for(t in seq_len(length(int.lolli$x$data[[n]]$text))){
            lipid <- NULL
            lipid <- stringr::str_split(
                int.lolli$x$data[[n]]$text[t], '<br />', n=2)[[1]][1] %>%
                stringr::str_remove('feature: ')
            int.lolli$x$data[[n]]$text[t] <- top.tab$hover[which(top.tab$feature == lipid)]
        }
    }
    return(list(static_de_lipid=lolli.plot,
                interactive_de_lipid=int.lolli,
                table_de_lipid=top.tab))
}


.dot_multiGroup <- function(deSp_result, lipid_char,
                            significant=c('pval', 'padj')){
    ## Color map
    class.color <- unique(c(RColorBrewer::brewer.pal(9, "Set1"),
                            RColorBrewer::brewer.pal(8, "Dark2"),
                            RColorBrewer::brewer.pal(8, "Set2"),
                            RColorBrewer::brewer.pal(8, "Accent"),
                            RColorBrewer::brewer.pal(12, "Set3"),
                            RColorBrewer::brewer.pal(9, "Pastel1"),
                            RColorBrewer::brewer.pal(8, "Pastel2")))

    ## Join deSp_result and lipid_char
    deSp.res <-  deSp_result %>% dplyr::left_join(lipid_char, by='feature') %>%
        dplyr::arrange(class, feature) %>%
        dplyr::mutate(index=dplyr::row_number())

    ## Dot plot
    if(significant == 'pval'){
        # Table
        dot.tab <- deSp.res %>% dplyr::mutate(
            hover=paste0(
                'Lipid: ', feature, '<br />Class: ', class,
                '<br />-log10(pval): ', format(negLog10pval, digits=3),
                '<br />pval: ', format(pval, digits=3, scientific=TRUE),
                '<br />Serial number: ', index),
            alpha=ifelse(sig_pval == 'yes', 1, 0.2))
        # Plot
        dot.plot <- ggpubr::ggscatter(
            data=dot.tab, x='index', y='negLog10pval', color='class', fill='class',
            alpha='alpha', palette=class.color, xlab='Lipid',
            ylab='-log10(pval)', legend='right') +
            ggplot2::guides(alpha='none')
    }else if(significant == 'padj'){
        # Table
        dot.tab <- deSp.res %>% dplyr::mutate(
            hover=paste0(
                'Lipid: ', feature, '<br />Class: ', class,
                '<br />-log10(padj): ', format(negLog10padj, digits=3),
                '<br />padj: ', format(padj, digits=3, scientific=TRUE),
                '<br />Serial number: ', index),
            alpha=ifelse(sig_padj == 'yes', 1, 0.2))
        # Plot
        dot.plot <- ggpubr::ggscatter(
            data=dot.tab, x='index', y='negLog10padj', color='class', fill='class',
            alpha='alpha', palette=class.color, xlab='Lipid',
            ylab='-log10(padj)', legend='right') +
            ggplot2::guides(alpha='none')
    }
    int.dot <- plotly::ggplotly(dot.plot)
    n <- t <- NULL
    for(n in seq_len(length(int.dot$x$data))){
        #int.dot$x$data[[n]]$marker$opacity <- 1
        for(t in seq_len(length(int.dot$x$data[[n]]$text))){
            idx <- stringr::str_split(
                int.dot$x$data[[n]]$text[t], '<br />', n=2)[[1]][1] %>%
                stringr::str_remove('.* ')
            int.dot$x$data[[n]]$text[t] <- dot.tab$hover[which(
                as.character(dot.tab$index) == idx)]
        }
    }
    return(list(static_dotPlot=dot.plot,
                interactive_dotPlot=int.dot,
                table_dotPlot=dot.tab))
}











