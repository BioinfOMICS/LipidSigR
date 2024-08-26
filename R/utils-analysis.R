## Two-groups analysis: Calculate mean, SD & FC
.mean_sd_twoGroup <- function(abundance, group_info){
    colnames(abundance)[1] <- 'feature'
    mean_sd_two <- abundance %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name') %>%
        dplyr::arrange(group, pair, feature) %>%
        dplyr::group_by(feature, group) %>%
        dplyr::summarise(value=list(abund), .groups='keep') %>%
        dplyr::mutate(mean=mean(unlist(value), na.rm=TRUE),
                      sd=sd(unlist(value), na.rm=TRUE)) %>%
        dplyr::mutate(mean=ifelse(is.nan(mean), NA, mean),
                      sd=ifelse(is.nan(sd), NA, sd)) %>%
        dplyr::select(-value) %>%
        tidyr::pivot_wider(names_from=group, values_from=c(mean, sd),
                           names_glue="{.value}_{group}") %>%
        dplyr::mutate(FC=mean_exp/mean_ctrl,
                      log2FC=log2(FC))
    return(mean_sd=mean_sd_two)
}

## Multi-groups analysis: Calculate mean & SD
.mean_sd_multiGroup <- function(abundance, group_info){
    colnames(abundance)[1] <- 'feature'
    mean_sd_multi <- abundance %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name') %>%
        dplyr::group_by(feature, group) %>%
        dplyr::summarise(value=list(abund), .groups='keep') %>%
        dplyr::mutate(mean=mean(unlist(value), na.rm=TRUE),
                      sd=sd(unlist(value), na.rm=TRUE)) %>%
        dplyr::mutate(mean=ifelse(is.nan(mean), NA, mean),
                      sd=ifelse(is.nan(sd), NA, sd)) %>%
        dplyr::select(-value) %>%
        tidyr::pivot_wider(names_from=group, values_from=c(mean, sd),
                           names_glue="{.value}_{group}")
    return(mean_sd=mean_sd_multi)
}


## Two-groups analysis: Statistical calculations
.stat_twoGroup <- function(abundance, group_info, paired_sample=FALSE,
                           test=c('t-test', 'Wilcoxon test')){
    colnames(abundance)[1] <- 'feature'
    ## gather the abundance data
    abund.ga <- abundance %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name') %>%
        dplyr::group_by(feature, group) %>%
        dplyr::summarise(value=list(abund), .groups='drop') %>%
        tidyr::spread(group, value) %>%
        dplyr::group_by(feature)
    ## t-test
    if(test == 't-test'){
        res.table <- abund.ga %>% dplyr::mutate(statistic=tryCatch(
            stats::t.test(unlist(exp), unlist(ctrl),
                          paired=paired_sample, var.equal=TRUE)$stat,
            error=function(e){NA}),
            pval=tryCatch(
                stats::t.test(unlist(exp), unlist(ctrl),
                              paired=paired_sample, var.equal=TRUE)$p.value,
                error=function(e){NA}))
    }else if(test == 'Wilcoxon test'){
        res.table <- abund.ga %>% dplyr::mutate(statistic=tryCatch(
            stats::wilcox.test(unlist(exp), unlist(ctrl),
                               paired=paired_sample)$stat,
            error=function(e){NA}),
            pval=tryCatch(
                stats::wilcox.test(
                    unlist(exp), unlist(ctrl), paired=paired_sample)$p.value,
                error=function(e){NA}))
    }
    res.table %<>% dplyr::mutate(method=test, .before='statistic') %>%
        dplyr::mutate(negLog10pval=-log10(pval)) %>%
        dplyr::select(-exp, -ctrl) %>% dplyr::ungroup()
    res.table$padj <- stats::p.adjust(res.table$pval, method='BH')
    res.table$negLog10padj <- -log10(res.table$padj)
    return(stat_res=res.table)
}


## Multi-groups analysis: Statistical calculations
.stat_multiGroup <- function(abundance, group_info,
                             test=c('One-way ANOVA', 'Kruskal–Wallis test')){
    colnames(abundance)[1] <- 'feature'
    # gather the abundance data
    abund.ga <- abundance %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(group_info, by='sample_name')

    res.table <- NULL
    for(i in 1:nrow(abundance)){
        sub.abund=test.res=test.table=NULL
        sub.abund <- abund.ga %>% dplyr::filter(feature == abundance$feature[i])
        ## One-way ANOVA
        if(test == 'One-way ANOVA'){
            test.res <- tryCatch({
                test.res <- stats::aov(abund ~ group, data=sub.abund)
            }, error = function(e) {
                test.res <- NULL
            })
            if(!is.null(test.res)) test.table <- broom::tidy(test.res)[1,]
        }else if(test == 'Kruskal–Wallis test'){ ## Kruskal–Wallis test
            test.res <- tryCatch({
                test.res <- stats::kruskal.test(abund ~ group, data=sub.abund)
            }, error = function(e) {
                test.res <- NULL
            })
            if(!is.null(test.res)) test.table <- broom::tidy(test.res)
        }
        ## statistic table
        if(!is.null(test.table)){
            if(!'statistic' %in% colnames(test.table)) test.table$statistic <- NA
            if(!'p.value' %in% colnames(test.table)) test.table$p.value <- NA
            test.table %<>%
                dplyr::mutate(
                    feature=abundance$feature[i],
                    negLog10pval=ifelse(!is.na(p.value), -log10(p.value), NA),
                    method=test) %>%
                dplyr::select(
                    feature, method, statistic, 'pval'='p.value', negLog10pval)
        }else{
            test.table <- data.frame(
                feature=abundance$feature[i], method=test, statistic=NA,
                pval=NA, negLog10pval=NA, stringsAsFactors=FALSE)
        }
        res.table <- data.table::rbindlist(list(res.table, test.table),
                                           use.names=TRUE, fill=TRUE)
    } #for loop
    res.table$padj <- stats::p.adjust(res.table$pval, method='BH')
    res.table$negLog10padj <- -log10(res.table$padj)
    return(stat_res=res.table)
}


## Two-way ANOVA
.two_way_anova <- function(abundance, group_info, char=NULL){
    ### abundance table
    ## first col: lipid characteristic
    ## second col: char feature
    colnames(abundance)[1:2] <- c('characteristic', 'feature')
    test.all <- NULL
    if(is.null(char)){
        charList <- unique(abundance$characteristic)
    }else{
        charList <- char
    }
    for(i in 1:length(charList)){
        CHAR=ABUND=test.res=test.table=NULL
        CHAR <- charList[i]
        ABUND <- abundance %>% dplyr::filter(characteristic == CHAR) %>%
            tidyr::gather(sample_name, abund, -1:-2) %>%
            dplyr::left_join(group_info, by='sample_name')
        test.res <- tryCatch({
            test.res <- stats::aov(abund ~ feature * group, data=ABUND)
        }, error = function(e){
            test.res <- NULL
        })
        if(!is.null(test.res)){
            test.table <- broom::tidy(test.res) %>%
                dplyr::mutate(characteristic=CHAR) %>%
                dplyr::filter(term != 'Residuals') %>%
                dplyr::select(characteristic, term, statistic, p.value) %>%
                tidyr::pivot_wider(
                    names_from='term', values_from=c('statistic', 'p.value'),
                    names_glue='{term}_{.value}') %>%
                dplyr::select(
                    characteristic, 'fval_2factors'='feature:group_statistic',
                    'pval_2factors'='feature:group_p.value',
                    'fval_feature'='feature_statistic',
                    'pval_feature'='feature_p.value',
                    'fval_group'='group_statistic',
                    'pval_group'='group_p.value')
        }else{
            test.table <- data.frame(characteristic=CHAR, fval_2factors=NA,
                                     pval_2factors=NA, fval_feature=NA,
                                     pval_feature=NA, fval_group=NA,
                                     pval_group=NA, stringsAsFactors=FALSE)
        }
        test.all <- data.table::rbindlist(list(test.all, test.table))
    }
    if(is.null(char)){
        test.all$padj_feature <- p.adjust(test.all$pval_feature, method='BH')
        test.all$padj_group <- p.adjust(test.all$pval_group, method='BH')
        test.all$padj_2factors <- p.adjust(test.all$pval_2factors, method='BH')
        test.all %<>% dplyr::select(
            characteristic, fval_2factors, pval_2factors, padj_2factors,
            fval_feature, pval_feature, padj_feature, fval_group, pval_group,
            padj_group)
    }
    return(twoWayAnova=test.all)
} #Function: .twoWayAnova


## Label and filter significant features
.sig_feature <- function(stat_table, significant=c('pval', 'padj', 'FC'),
                         p_cutoff=0.05, FC_cutoff=NULL){
    ## stat_table must have contain 'pval' and 'padj' columns,
    ## column 'log2FC' is optional.
    if(!is.null(p_cutoff) & !is.null(FC_cutoff)){
        all_table <- stat_table %>%
            dplyr::mutate(
                sig_pval=ifelse(pval < p_cutoff & abs(log2FC) > log2(FC_cutoff), 'yes', 'no'),
                sig_padj=ifelse(padj < p_cutoff & abs(log2FC) > log2(FC_cutoff), 'yes', 'no'))
    }else if(!is.null(p_cutoff) & is.null(FC_cutoff)){
        all_table <- stat_table %>%
            dplyr::mutate(
                sig_pval=ifelse(pval < p_cutoff, 'yes', 'no'),
                sig_padj=ifelse(padj < p_cutoff, 'yes', 'no'))
    }else if(is.null(p_cutoff) & !is.null(FC_cutoff)){
        all_table <- stat_table %>%
            dplyr::mutate(
                sig_FC=ifelse(abs(log2FC) > log2(FC_cutoff), 'yes', 'no'))
    }
    switch(significant,
           pval=sig_table <- all_table %>% dplyr::filter(sig_pval == 'yes') %>%
               dplyr::arrange(pval),
           padj=sig_table <- all_table %>% dplyr::filter(sig_padj == 'yes') %>%
               dplyr::arrange(padj),
           FC=sig_table <- all_table %>% dplyr::filter(sig_FC == 'yes') %>%
               dplyr::arrange(dplyr::desc(log2FC)))
    return(list(all_table=all_table,
                sig_table=sig_table))
}

## how to use
# if (.char_dataType(lipid_char_table, char)!="numeric") {
#    warning("Line plot and box plot will not produce for character lipid characteristic.")
# }
.char_dataType <- function(lipid_char_table, char) {
    char_type <- ifelse(isTRUE(class(lipid_char_table[, char]) %in% c("numeric", "integer")),
           "numeric", "char")
    return(char_type=char_type)
}

## Network species to class
.nw_classInfo <- function(abundance, lipid_char){
    ## Reference network species node
    spNode <- networkNode %>% dplyr::filter(!is.na(sum_composition), !is.na(goslin))
    ## Lipid class-network ID pair
    class.networkID.pair <- lipid_char %>%
        dplyr::distinct(class, LIPIDMAPS.reaction.abbr)
    ## Summarize abundance by lipid class
    summ.abund <- abundance %>%
        tidyr::gather(sample_name, abund, -1) %>%
        dplyr::left_join(lipid_char[,c('feature', 'class', 'LIPIDMAPS.reaction.abbr')], by='feature') %>%
        dplyr::group_by(LIPIDMAPS.reaction.abbr, sample_name) %>%
        dplyr::summarise(abund=sum(abund, na.rm=TRUE), .groups='drop') %>%
        dplyr::filter(!is.na(LIPIDMAPS.reaction.abbr)) %>%
        dplyr::left_join(class.networkID.pair, by='LIPIDMAPS.reaction.abbr') %>%
        dplyr::left_join(spNode[,c('id', 'label')], by=c('LIPIDMAPS.reaction.abbr'='id')) %>%
        tidyr::spread(sample_name, abund) %>%
        dplyr::mutate(class=ifelse(!is.na(label), label, class)) %>%
        dplyr::select(-label)
    ## Class abundance
    class.abund <- summ.abund %>%
        dplyr::select('feature'='class', dplyr::everything(), -LIPIDMAPS.reaction.abbr)
    ## Mapped to network ID
    network.id <- summ.abund %>%
        dplyr::select('feature'='class', 'networkId'='LIPIDMAPS.reaction.abbr')
    return(list(class.abund=class.abund,
                network.id=network.id))
}

#### Calculate characteristic index of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.char_index <- function(abundance, char){
    abund <- abundance %>% dplyr::mutate(feature=as.numeric(feature))
    combined.abund <- purrr::map2(
        abund[1], abund[-1], ~sum(.x*.y, na.rm = T)/sum(.y, na.rm = T)) %>%
        as.data.frame() %>% dplyr::mutate(index=paste0(char, ' index')) %>%
        dplyr::select(index, everything())
    colnames(combined.abund) <- colnames(abund)
    return(combined.abund)
}

#### Multi-groups analysis: plot table of bar plot and line plot of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.deChar_plot_tab_multiGroup <- function(stat_tab, group_info, post_hoc_sig){
    if(post_hoc_sig == 'pval'){
        stat_tab$p <- stat_tab$post_hoc_pval
        stat_tab$significance <- stat_tab$post_hoc_sig_pval
    }else if(post_hoc_sig == 'padj'){
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
                '<br />', post_hoc_sig,': ', format(p, digits=3, scientific=TRUE),
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

#### Multi-groups analysis: Static bar plot of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.deChar_bar_multiGroup <- function(plot_tab, char, star_position, plot_title=NULL){
    # Bar plot
    bar <- ggplot2::ggplot(
        data=plot_tab, ggplot2::aes(x=feature, y=mean, fill=group)) +
        ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
        ggplot2::scale_fill_brewer(palette = 'Set2') +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin=mean, ymax=mean+sd), color="gray39", width=.9,
            position=ggplot2::position_dodge()) +
        ggplot2::geom_text(ggplot2::aes(x=feature, y=max_error_bar+star_position,
                                        label=pvalue_text), color="red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=char, y='Abundance mean', title=plot_title) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=16, hjust=0.5),
                       axis.text=ggplot2::element_text(size=14),
                       axis.title=ggplot2::element_text(size=16),
                       legend.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=14))
    bar.plotly <- .deChar_bar_multiGroup_plotly(plot_tab, bar)
    # Sqrt bar plot
    bar.sqrt <- bar + ggplot2::scale_y_sqrt()
    bar.sqrt.plotly <- .deChar_bar_multiGroup_plotly(plot_tab, bar.sqrt)
    return(list(bar.tab=plot_tab,
                bar.plot=bar,
                bar.plotly=bar.plotly,
                bar.sqrt.plot=bar.sqrt,
                bar.sqrt.plotly=bar.sqrt.plotly))
}

#### Multi-groups analysis: Interactive bar plot of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.deChar_bar_multiGroup_plotly <- function(plot_tab, bar){
    bar.plotly <- plotly::ggplotly(bar)
    for(i in 1:length(bar.plotly$x$data)){
        for(j in 1:length(bar.plotly$x$data[[i]]$text)){
            if(i %in% 1:(length(unique(plot_tab$group))*2)){
                legendgroup <- bar.plotly$x$data[[i]]$legendgroup
                f <- stringr::str_split(
                    bar.plotly$x$data[[i]]$text[j], pattern='<br />')[[1]][1] %>%
                    stringr::str_remove('feature: ')
                bar.plotly$x$data[[i]]$text[j] <- plot_tab$hover[which(plot_tab$feature == f
                                                                       & plot_tab$group == legendgroup)]
            }else{
                f <- stringr::str_split(
                    bar.plotly$x$data[[i]]$hovertext[j], pattern='<br />')[[1]][1] %>%
                    stringr::str_remove('feature: ')
                pt <- stringr::str_split(
                    bar.plotly$x$data[[i]]$hovertext[j], pattern='<br />')[[1]][3] %>%
                    stringr::str_remove('pvalue_text: ')
                bar.plotly$x$data[[i]]$hovertext[j] <- plot_tab$hovertext[which(plot_tab$feature == f
                                                                                & plot_tab$pvalue_text == pt)][1]
            }
        }
    }
    return(bar.plotly)
}

#### Multi-groups analysis: Static line plot of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.deChar_line_multiGroup <- function(plot_tab, char, star_position, plot_title=NULL){
    # Line plot
    line <- ggplot2::ggplot(
        data=plot_tab, ggplot2::aes(x=feature, y=mean, group=group, color=group)) +
        ggplot2::geom_line(stat="identity", position=ggplot2::position_dodge(0.05))+
        ggplot2::scale_color_brewer(palette = 'Set2') +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin=mean, ymax=mean+sd), color="gray39",
            position=ggplot2::position_dodge(0.05)) +
        ggplot2::geom_text(
            ggplot2::aes(x=feature, y=max_error_bar+star_position,
                         label=pvalue_text), color="red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x=char, y='Abundance mean', title=plot_title) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=16, hjust=0.5),
                       axis.text=ggplot2::element_text(size=14),
                       axis.title=ggplot2::element_text(size=16),
                       legend.text=ggplot2::element_text(size=14),
                       legend.title=ggplot2::element_text(size=14))
    line.plotly <- .deChar_line_multiGroup_plotly(plot_tab, line)
    # Sqrt bar plot
    line.sqrt <- line + ggplot2::scale_y_sqrt()
    line.sqrt.plotly <- .deChar_line_multiGroup_plotly(plot_tab, line.sqrt)
    return(list(line.tab=plot_tab,
                line.plot=line,
                line.plotly=line.plotly,
                line.sqrt.plot=line.sqrt,
                line.sqrt.plotly=line.sqrt.plotly))
}

#### Multi-groups analysis: Interactive line plot of DE characteristics ####
## Used functions: plot_deChar_multiGroup(), plot_subChar_multiGroup()
.deChar_line_multiGroup_plotly <- function(plot_tab, line.plot){
    line.plotly <- plotly::ggplotly(line.plot)
    for(i in 1:length(line.plotly$x$data)){
        if(i == (length(unique(plot_tab$group))+1)){
            line.plotly$x$data[[i]]$text <- ""
        }else{
            for(j in 1:length(line.plotly$x$data[[i]]$text)){
                if(i %in% 1:length(unique(plot_tab$group))){
                    legendgroup <- line.plotly$x$data[[i]]$legendgroup
                    f <- stringr::str_split(
                        line.plotly$x$data[[i]]$text[j], pattern='<br />')[[1]][1] %>%
                        stringr::str_remove('feature: ')
                    line.plotly$x$data[[i]]$text[j] <- plot_tab$hover[which(plot_tab$feature == f
                                                                            & plot_tab$group == legendgroup)]
                }else{
                    f <- stringr::str_split(
                        line.plotly$x$data[[i]]$hovertext[j], pattern='<br />')[[1]][1] %>%
                        stringr::str_remove('feature: ')
                    pt <- stringr::str_split(
                        line.plotly$x$data[[i]]$hovertext[j], pattern='<br />')[[1]][3] %>%
                        stringr::str_remove('pvalue_text: ')
                    line.plotly$x$data[[i]]$hovertext[j] <- plot_tab$hovertext[which(plot_tab$feature == f
                                                                                     & plot_tab$pvalue_text == pt)][1]
                }
            }
        }
    }
    return(line.plotly)
}

#### Multi-groups analysis: statistical test & post hoc test for plotting box plot. ####
## Used functions: boxPlot_feature_multiGroup(), plot_deChar_multiGroup(), heatmap_chain_db()
.box_stat_multiGroup <- function(sub_abund, group_info, ref_group, test){
    test.res=test.table=post.res=post.table=NULL
    post.hoc.method <- ifelse(test == 'One-way ANOVA', "Tukey's HSD", "Dunn's Test")
    if(test == 'One-way ANOVA'){
        test.res <- tryCatch({
            test.res <- stats::aov(abund ~ original_group_name, data=sub_abund)
        }, error=function(e) {
            test.res <- NULL
        })
        if(!is.null(test.res)){
            test.table <- broom::tidy(test.res)[1,]
            if(!'statistic' %in% colnames(test.table)) test.table$statistic <- NA
            if(!'p.value' %in% colnames(test.table)) test.table$p.value <- NA
            if(test.table$p.value < 0.05){
                post.res <- tryCatch({
                    post.res <- rstatix::tukey_hsd(
                        sub_abund, abund ~ original_group_name, detailed=TRUE)
                }, warning=function(w){
                    post.res <- NULL
                }, error=function(e) {
                    post.res <- NULL
                })
            }
        }
    }else if(test == 'Kruskal–Wallis test'){
        test.res <- tryCatch({
            test.res <- stats::kruskal.test(abund ~ original_group_name,
                                            data=sub_abund)
        }, error=function(e){
            test.res <- NULL
        })
        if(!is.null(test.res)){
            test.table <- broom::tidy(test.res)
            if(!'statistic' %in% colnames(test.table)) test.table$statistic <- NA
            if(!'p.value' %in% colnames(test.table)) test.table$p.value <- NA
            if(test.table$p.value < 0.05){
                post.res <- tryCatch({
                    post.res <- rstatix::dunn_test(
                        data=sub_abund, abund ~ original_group_name,
                        p.adjust.method='BH')
                }, warning=function(w){
                    post.res <- NULL
                }, error=function(e){
                    post.res <- NULL
                })
            }
        }
    }

    if(is.null(test.res)){
        comparison <- combn(sort(unique(sub_abund$original_group_name)), m=2, simplify=FALSE)
        com.all <- NULL
        for(i in 1:length(comparison)){
            com <- paste0(comparison[[i]], collapse='-')
            com.all <- c(com.all, com)
        }
        res.table <- data.frame(feature=unique(sub_abund$feature),
                                contrast=com.all, method=test, statistic=NA,
                                pval=NA, post_hoc_method=NA,
                                post_hoc_pval=NA, post_hoc_padj=NA,
                                stringsAsFactors=FALSE)
        if(post.hoc.method == "Tukey's HSD"){
            res.table %<>% dplyr::select(-post_hoc_pval)
        }
        return(list(res_table=res.table,
                    post_hoc_table=NULL))
    }

    if(!is.null(post.res)){
        if(post.hoc.method == "Tukey's HSD") post.res %<>% dplyr::mutate(p=NA)
        post.table <- post.res %>%
            dplyr::mutate(contrast=paste0(group1, '-', group2),
                          post_hoc_method=post.hoc.method) %>%
            dplyr::select(contrast, post_hoc_method, 'post_hoc_pval'='p',
                          'post_hoc_padj'='p.adj')
        test.table$feature <- unique(sub_abund$feature)
        test.table$method <- test
        test.table %<>% dplyr::select(feature, statistic, 'pval'='p.value', method)
        res.table <- cbind(test.table, post.table) %>%
            dplyr::select(feature, contrast, method, statistic, pval,
                          post_hoc_method, post_hoc_pval, post_hoc_padj) %>%
            dplyr::arrange(post_hoc_pval)
        if(post.hoc.method == "Tukey's HSD"){
            res.table %<>% dplyr::select(-post_hoc_pval)
            post.res %<>% dplyr::select(-p)
        }
    }else{
        comparison <- combn(sort(unique(sub_abund$original_group_name)),
                            m=2, simplify=FALSE)
        com.all <- NULL
        for(i in 1:length(comparison)){
            com <- paste0(comparison[[i]], collapse='-')
            com.all <- c(com.all, com)
        }
        res.table <- data.frame(
            feature=unique(sub_abund$feature), contrast=com.all, method=test,
            statistic=test.table$statistic, pval=test.table$p.value,
            post_hoc_method=NA, post_hoc_pval=NA, post_hoc_padj=NA,
            stringsAsFactors=FALSE)
    }

    return(list(res_table=res.table,
                post_hoc_table=post.res))
}

#### Multi-groups analysis: static box plot with p-values. ####
## Used functions: boxPlot_feature_multiGroup(), plot_deChar_multiGroup(), heatmap_chain_db()
.box_plot_multiGroup <- function(
        sub_abund, group_info, test, post_hoc_sig, res_table,
        post_hoc_table, plot_title, plot_ylabel){

    post.hoc.method <- ifelse(test == 'One-way ANOVA', "Tukey's HSD", "Dunn's Test")
    fac.level <- group_info %>%
        dplyr::distinct(original_group_name) %>%
        dplyr::mutate(x=1:nrow(.))
    if(!is.null(post_hoc_table)){
        # Autocompute P-value Positions For Plotting Significance
        pwc <- post_hoc_table %>%
            rstatix::add_xy_position(x='original_group_name') %>%
            dplyr::left_join(fac.level, by=c('group1'='original_group_name')) %>%
            dplyr::left_join(fac.level, by=c('group2'='original_group_name')) %>%
            dplyr::mutate(xmin=x.x,
                          xmax=x.y) %>%
            dplyr::select(-x.x, -x.y)
        # Add p-value significance symbols into a data frame
        if(post.hoc.method == "Dunn's Test" && post_hoc_sig == 'pval'){
            pwc %<>% rstatix::add_significance(p.col='p') %>%
                dplyr::mutate(p.adj.signif=p.signif)
        }
    }
    # Group as factor
    sub_abund$original_group_name <- factor(sub_abund$original_group_name,
                                            levels=fac.level$original_group_name)
    # Box plot
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

    if(!is.null(post_hoc_table)){
        BOX <- BOX +
            ggpubr::stat_pvalue_manual(
                pwc, label='p.adj.signif', hide.ns=TRUE, size=6, bracket.size=0.5) +
            ggplot2::labs(caption=rstatix::get_pwc_label(pwc))
    }
    return(list(box.plot=BOX,
                box.tab=sub_abund))
}
