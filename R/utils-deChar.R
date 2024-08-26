.table_deChar_twoGroup <- function(stat_tab, group_info, significant){
    if(significant == 'pval'){
        stat_tab$p <- stat_tab$post_hoc_pval
        stat_tab$significance <- stat_tab$post_hoc_sig_pval
    }else if(significant == 'padj'){
        stat_tab$p <- stat_tab$post_hoc_padj
        stat_tab$significance <- stat_tab$post_hoc_sig_padj
    }
    colName <- c(colnames(stat_tab)[stringr::str_starts(colnames(stat_tab), 'mean_')],
                 colnames(stat_tab)[stringr::str_starts(colnames(stat_tab), 'sd_')])
    plot_tab <- stat_tab %>%
        dplyr::select(characteristic, feature, significance, p, all_of(colName)) %>%
        tidyr::pivot_longer(
            cols=dplyr::all_of(colName), cols_vary='slowest',
            names_to=c('.value', 'group'), names_pattern='(.*)_(.*)') %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(max_error_bar=max(mean+sd)) %>% dplyr::ungroup() %>%
        dplyr::mutate(
            pvalue_text=ifelse(p < 0.001 ,"***", ifelse(p < 0.01 ,"**", ifelse(p < 0.05 ,"*",""))),
            hover=paste0(
                'Characteristic: ', characteristic, '<br />Feature: ', feature,
                '<br />Group: ', group, '<br />Mean: ', format(mean, digits=3),
                '<br />SD: ', format(sd, digits=3)),
            hovertext=paste0(
                'Characteristic: ', characteristic, '<br />Feature: ', feature,
                '<br />', significant,': ', format(p, digits=3, scientific=TRUE),
                '<br />Significance: ', significance))
    if(unique(plot_tab$characteristic) %in%
       c('Total.C', 'Total.DB', 'Total.OH', 'FA.C', 'FA.DB', 'FA.OH')){
        plot_tab$feature <- as.factor(as.numeric(plot_tab$feature))
    }else{
        plot_tab$feature <- as.factor(plot_tab$feature)
    }
    fac.level <- group_info %>% dplyr::distinct(group)
    plot_tab$group <- factor(plot_tab$group, levels = fac.level$group)
    return(plot_tab)
}


.barPlot_deChar_twoGroup <- function(plot_tab, char, split_class=NULL) {
    star_number <- ifelse(is.null(split_class), 5, 3)
    barPlot <- ggplot2::ggplot(
        data=plot_tab,
        ggplot2::aes(x=feature, y=mean, fill=group)) +
        ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
        ggplot2::scale_fill_manual(values=c('lightslateblue', 'sienna2')) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin=mean, ymax=mean+sd), color="gray39", width=.9,
            position=ggplot2::position_dodge()) +
        ggplot2::geom_text(ggplot2::aes(
                x=feature, y=max_error_bar+star_number, label=pvalue_text), color="red") + ## bar
        ggplot2::theme_minimal() +
        ggplot2::labs(x=char, title=split_class) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=14, hjust = 0.5))
    in_barPlot <- plotly::ggplotly(barPlot)
    for (i in seq_len(length(unique(barPlot$data$group)))) {
        n <- length(unique(barPlot$data$group))
        data <- barPlot$data[which(barPlot$data$group == unique(barPlot$data$group)[i]),]
        in_barPlot$x$data[[i]]$text <- paste0(
            "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        in_barPlot$x$data[[i+n]]$text <- paste0(
            "feature : ", data$feature, "\nmean :", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
    }
    for (i in seq_len(length(in_barPlot$x$data))) {
        text <- stringr::str_split(in_barPlot$x$data[[i]]$hovertext, "<br />max_error_bar")
        hovertext <- list()
        if(length(text) > 0){
            for (j in seq_len(length(text))) {
                hovertext[[j]] <- paste(text[[j]][1], "\nsignificant : YES")
            }
            in_barPlot$x$data[[i]]$hovertext <- hovertext
        }
    }
    ## sqrt
    barPlot_sqrt <- barPlot + ggplot2::scale_y_sqrt()
    in_barPlot_sqrt <- plotly::ggplotly(barPlot_sqrt)
    for (i in seq_len(length(unique(barPlot_sqrt$data$group)))) {
        n <- length(unique(barPlot_sqrt$data$group))
        data <- barPlot_sqrt$data[which(barPlot_sqrt$data$group == unique(barPlot_sqrt$data$group)[i]),]
        in_barPlot_sqrt$x$data[[i]]$text <- paste0(
            "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        in_barPlot_sqrt$x$data[[i+n]]$text <- paste0(
            "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
    }
    for (i in seq_len(length(in_barPlot_sqrt$x$data))) {
        text <- stringr::str_split(in_barPlot_sqrt$x$data[[i]]$hovertext, "<br />max_error_bar")
        hovertext <- list()
        if(length(text) > 0){
            for (j in seq_len(length(text))) {
                hovertext[[j]] <- paste(text[[j]][1], "\nsignificant : YES")
            }
            in_barPlot_sqrt$x$data[[i]]$hovertext <- hovertext
        }
    }
    return(list(
        in_barPlot=in_barPlot, barPlot=barPlot,
        in_barPlot_sqrt=in_barPlot_sqrt, barPlot_sqrt=barPlot_sqrt))
}


.lintPlot_deChar_twoGroup <- function(plot_tab, char, split_class=NULL) {
    star_number <- ifelse(is.null(split_class), 5, 3)

    linePlot <- ggplot2::ggplot(
        data=plot_tab, ggplot2::aes(x=feature, y=mean, group=group, color=group)) +
        ggplot2::geom_line(stat="identity", position=ggplot2::position_dodge(0.05))+
        ggplot2::scale_color_manual(values=c('lightslateblue', 'sienna2')) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin=mean, ymax=mean+sd), color="gray39",
            position=ggplot2::position_dodge(0.05)) +
        ggplot2::geom_text(
            ggplot2::aes(x=feature, y=max_error_bar+star_number,label=pvalue_text), color="red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=char, title=split_class) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=14, hjust = 0.5))
    in_linePlot <- plotly::ggplotly(linePlot)
    if (!is.null(split_class)) {
        for(i in seq_len(length(in_linePlot$x$data))){
            if(!is.null(in_linePlot$x$data[i]$name)){
                in_linePlot$x$data[i]$name <- gsub("\\(", "", stringr::str_split(in_linePlot$x$data[i]$name, ",")[[1]][1])
            }
        }
    }
    for (i in seq_len(length(unique(linePlot$data$group)))) {
        n <- length(unique(linePlot$data$group))
        data <- linePlot$data[which(linePlot$data$group == unique(linePlot$data$group)[i]),]
        in_linePlot$x$data[[i]]$text <- paste0(
            "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        if(sum(grepl("\\*",in_linePlot$x$data[[i+n]]$text)) == 0){
            in_linePlot$x$data[[i+n]]$text <- paste0(
                "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
                "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        }
    }
    for (i in seq_len(length(in_linePlot$x$data))) {
        text <- stringr::str_split(in_linePlot$x$data[[i]]$hovertext, "<br />max_error_bar")
        hovertext <- list()
        if(length(text) > 0){
            for (j in seq_len(length(text))) {
                hovertext[[j]] <- paste(text[[j]][1], "\nsignificant : YES")
            }
            in_linePlot$x$data[[i]]$hovertext <- hovertext
        }
    }
    ## sqrt
    linePlot_sqrt <- linePlot + ggplot2::scale_y_sqrt()
    in_linePlot_sqrt <- plotly::ggplotly(linePlot_sqrt)
    if (!is.null(split_class)) {
        for(i in seq_len(length(in_linePlot_sqrt$x$data))){
            if(!is.null(in_linePlot_sqrt$x$data[i]$name)){
                in_linePlot_sqrt$x$data[i]$name <- gsub("\\(", "", stringr::str_split(in_linePlot_sqrt$x$data[i]$name, ",")[[1]][1])
            }
        }
    }
    for (i in seq_len(length(unique(linePlot$data$group)))) {
        n <- length(unique(linePlot$data$group))
        data <- linePlot$data[which(linePlot$data$group == unique(linePlot$data$group)[i]),]
        in_linePlot_sqrt$x$data[[i]]$text <- paste0(
            "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
            "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        if(sum(grepl("\\*",in_linePlot_sqrt$x$data[[i+n]]$text)) == 0){
            in_linePlot_sqrt$x$data[[i+n]]$text <- paste0(
                "feature : ", data$feature, "\nmean : ", round(data$mean, 3),
                "\nsd : ", round(data$sd, 3), "\ngroup : ", data$group)
        }
    }
    for (i in seq_len(length(in_linePlot_sqrt$x$data))) {
        text <- stringr::str_split(in_linePlot_sqrt$x$data[[i]]$hovertext, "<br />max_error_bar")
        hovertext <- list()
        if(length(text) > 0){
            for (j in seq_len(length(text))) {
                hovertext[[j]] <- paste(text[[j]][1], "\nsignificant : YES")
            }
            in_linePlot_sqrt$x$data[[i]]$hovertext <- hovertext
        }
    }
    return(list(linePlot=linePlot, in_linePlot=in_linePlot,
                linePlot_sqrt=linePlot_sqrt, in_linePlot_sqrt=in_linePlot_sqrt))
}

.boxPlot_deChar_twoGroup <- function(Combine_char_result_table, boxTab, char, split_class=NULL) {
    if(isTRUE(Combine_char_result_table$sig == 'yes') ){
        group_name <- c(
            unique(boxTab$group)[1], paste0(unique(boxTab$group)[1], "0"),
            unique(boxTab$group)[2])
        group_max <- max(boxTab$abund) %>% unique()
        group_min <- min(boxTab$abund) %>% unique()
        if (Combine_char_result_table[[significant]] <= 0.05 & Combine_char_result_table[[significant]] > 0.01){
            t.text <- c('', '*', '')
        }else if (Combine_char_result_table[[significant]] <= 0.01 & Combine_char_result_table[[significant]] > 0.001){
            t.text <- c('', '**', '')
        }else{
            t.text <- c('', '***','')
        }
        boxTab_ctrl <- boxTab %>% dplyr::filter(group == "ctrl")
        boxTab_exp <- boxTab %>% dplyr::filter(group == "exp")
        in_boxPlot <- plotly::plot_ly() %>%
            plotly::add_bars(
                x=group_name, y=rep(group_max+0.15, 3), opacity=1, showlegend=FALSE,
                marker=list(line=list(color='rgba(0,0,0,0)'), color='rgba(0,0,0,0)'),
                textfont=list(color='red'), text=t.text, hoverinfo='none',
                textposition='outside', legendgroup="1") %>%
            plotly::add_lines(
                x=c(rep(paste("ctrl"), 2),
                    rep(paste("exp"), 2)),
                y=c(group_max+0.1, group_max+0.15, group_max+0.15, group_max+0.1),
                showlegend=FALSE, line=list(color='black'), legendgroup="1",
                hoverinfo='none') %>%
            plotly::add_boxplot(
                data=boxTab_ctrl, x=~group, y=~abund, color=I('lightslateblue'),
                name="ctrl", boxpoints="all", jitter=0.85,
                pointpos=0, marker=list(size=5, opacity=0.8)) %>%
            plotly::add_boxplot(
                data=boxTab_exp, x=~group, y=~abund, color=I('sienna2'),
                name="exp", boxpoints="all", jitter=0.85,
                pointpos=0, marker=list(size=5, opacity=0.8)) %>%
            plotly::layout(
                title=split_class, xaxis=list(
                    title='Group', tickmode='array',
                    tickvals= c("ctrl", '', "exp"),
                    ticktext=c("ctrl", '', "exp"),
                    titlefont=list(size=16), tickfont=list(size=14)),
                yaxis=list(title=paste0(char, ' index'), titlefont=list(size=16),
                           tickfont=list(size=14), range=c(group_min, group_max+0.5)),
                legend=list(font=list(size=14), y=0.5),
                margin=list(l=70, r=70, b=80, t=60))
        boxPlot <- ggpubr::ggboxplot(
            boxTab, x = "group", y = "abund", color = "group", add = "jitter") +
            ggplot2::scale_color_manual(values=c("lightslateblue", "sienna2")) +
            ggpubr::stat_compare_means(method = "t.test") +
            ggplot2::labs(
                y=paste0(char, ' index'), x='Group', title=split_class) +
            ggplot2::guides(color="none")
    }else{
        ## pvalNA
        in_boxPlot <- plotly::plot_ly(
            data=boxTab, x=~group, y=~abund, type='box', color=~group,
            colors=c('lightslateblue', 'sienna2'), boxpoints='all', jitter=0.85,
            pointpos=0, marker=list(size=5, opacity=0.8)) %>%
            plotly::layout(title=split_class, xaxis=list(
                title='Group', titlefont=list(size=16), tickfont=list(size=14)),
                yaxis=list(
                    title=paste0(char, ' index'), titlefont=list(size=16),
                    tickfont=list(size=14)),
                legend=list(font=list(size=14), y=0.5),
                margin=list(l=70, r=70, b=80, t=60))
        boxPlot <- ggpubr::ggboxplot(
            boxTab, x="group", y="abund", color="group", add="jitter") +
            ggplot2::scale_color_manual(values=c("lightslateblue", "sienna2")) +
            ggplot2::labs(
                y=paste0(char, ' index'), x='Group', title=split_class) +
            ggplot2::guides(color="none")
    }
    return(list(in_boxPlot=in_boxPlot, boxPlot=boxPlot))
}
