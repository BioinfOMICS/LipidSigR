#' @title deSp_twoGroup
#' @description Compute differentially expressed analysis of two groups
#' (independent) to find significant lipid species.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param ref_group Character. Group name of the reference group. It must be one of the
#' group names in the group information table's group column.
#' @param test Character. The method to use for comparing means. Allowed method include "t-test" and
#' "Wilcoxon test". Default is \code{'t-test'}.
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
#'     normalization='Percentage')
#' deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
#'     significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

deSp_twoGroup <- function(
        processed_se, ref_group, test=c('t-test', 'Wilcoxon test'),
        significant=c('pval', 'padj'),
        p_cutoff=0.05, FC_cutoff=1, transform=c('none', 'log10', 'square', 'cube')){
    ## check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if (is.null(test) | isFALSE(test %in% c('t-test', 'Wilcoxon test')) ) {
        stop("test must be one of 't-test' or 'Wilcoxon test''.")
    }
    if (is.null(significant) | isFALSE(significant %in% c('pval', 'padj')) ) {
        stop("significant must be one of 'pval' or 'padj'.")
    }
    if (!is.numeric(p_cutoff) | isFALSE(.check_numeric_range(p_cutoff, 0, 1)) ) {
        stop("p_cutoff must be a numeric value between 0 and 1.")
    }
    if (!is.numeric(FC_cutoff) | isFALSE(.check_numeric_range(FC_cutoff, 1, NULL)) ) {
        stop("FC_cutoff must be a numeric value >= 1.")
    }
    if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
        stop("transform must be one of 'none', 'log10', 'square', or 'cube'.")
    }

    abundance <- .extract_df(processed_se, type = "abundance")
    processed_abund <- abundance
    group_info_raw <- .extract_df(processed_se, type = "group")
    lipid_char_table <- .extract_df(processed_se, type = "lipid")
    .check_imputation(abundance)
    if (isTRUE(.check_nor_negative(abundance[-1])) ) {
        stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }
    if (.check_nGroup(group_info_raw) != 'two') {
        stop('This function is used for two groups. If your data consists of multi groups, please use deSp_multiGroup function for analysis.')
    }
    if (is.null(ref_group) | isFALSE(ref_group %in% unique(group_info_raw$group)) ){
        stop("ref_group must be one of the group names in the 'group' column of the group information table.")
    }
    ## add original group names
    group_info <- group_info_raw %>%
        dplyr::mutate("original_group_name"=group_info_raw$group, .before=group)
    group_info[group_info$original_group_name==ref_group, "group"] <- "ctrl"
    group_info[group_info$original_group_name!=ref_group, "group"] <- "exp"
    ## paired sample
    paired_sample <- ifelse(isFALSE(all(is.na(group_info$pair))), TRUE, FALSE)
    ## mean_ctrl, mean_exp, method, FC, log2FC
    abundance_tab <- .mean_sd_twoGroup(abundance, group_info)
    ## Data transformation
    abundance <- .transform(abundance, transform)
    if (sum(group_info$group=="ctrl")==1 && sum(group_info$group=="exp")==1) {
        res_table <- abundance_tab %>% dplyr::select(-c("sd_ctrl", "sd_exp"))
        diff_table <- .sig_feature(res_table, significant='FC', p_cutoff=NULL, FC_cutoff=FC_cutoff)
    } else {
        ## statistic table
        stat_table <- .stat_twoGroup(abundance, group_info, paired_sample, test)
        ## result table
        res_table <- abundance_tab %>% dplyr::left_join(stat_table, by='feature') %>%
            as.data.frame()
        diff_table <- .sig_feature(res_table, significant, p_cutoff, FC_cutoff)
    }
    diff_table_sig <- diff_table$sig_table
    if (nrow(diff_table_sig)==0) {
        warning("This case does not include any significant lipids, which will prevent some subsequent analyses from being performed.")
        sig_table <- "No significant lipids."
    } else if (sum(group_info$group=="ctrl")==1 && sum(group_info$group=="exp")==1) {
        sig_table <- diff_table_sig %>% dplyr::select(feature, sig_FC, everything())
    } else {
        sig_table <- diff_table_sig %>%
            dplyr::select(feature, paste0("sig_", significant), everything())
    }
    #exp_stat <- abundance_stat %>% dplyr::left_join(lipid_char_table, by="feature")

    abundance_mat <- abundance %>% dplyr::arrange(feature) %>%
        tibble::column_to_rownames(var="feature")
    lipid_char_table %<>% dplyr::arrange(feature)
    deSp_se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(abundance=as.matrix(abundance_mat)),
        rowData=S4Vectors::DataFrame(lipid_char_table),
        colData=S4Vectors::DataFrame(group_info),
        metadata=list(all_deSp_result=diff_table$all_table,
                      sig_deSp_result=sig_table,
                      processed_abundance=processed_abund,
                      significant=significant,
                      p_cutoff=p_cutoff,
                      FC_cutoff=FC_cutoff,
                      transform=transform))
    ## Check results
    # if(nrow(res.tab$sig_table) == 0){
    #     warning('')
    # }

    ## ctrl / exp
    # original_ctrl <- unique(group_info[group_info$group=="ctrl", "original_group_name"])
    # original_exp <- unique(group_info[group_info$group=="exp", "original_group_name"])
    # all_table_web <- abundance_stat %>% dplyr::mutate(
    #     comparison = paste0(original_exp, "(exp) / ", original_ctrl, "(ctrl)"))
    return(deSp_se)
}

#' @title plot_deSp_twoGroup
#' @description This function is for plotting the results of lipid species differential expression analysis.
#' @param deSp_se A SummarizedExperiment object with results computed by \code{\link{deSp_twoGroup}}.
#' @return Return a list of 3 interactive plots, 3 static plots, and 2 tables.
#' \enumerate{
#' \item interactive_de_lipid & static_de_lipid: a lollipop chart reveals the lipid species that pass chosen cut-off.
#' \item interactive_maPlot & static_maPlot: an MA plot of lipid species.
#' \item interactive_volcanoPlot & static_volcanoPlot: a volcano plot of lipid species.
#' \item table_de_lipid: table for plotting DE lipid plot.
#' \item table_ma_volcano: table for plotting MA plot and volcano plot.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
#'     significant='padj', p_cutoff=0.05, FC_cutoff=1, transform='log10')
#' deSp_plot <- plot_deSp_twoGroup(deSp_se)
#'
plot_deSp_twoGroup <- function(deSp_se){
    .check_de_outputSE(deSp_se, de_type="deSp")
    group_info <- .extract_df(deSp_se, type="group")
    all_table <- S4Vectors::metadata(deSp_se)[["all_deSp_result"]]
    sig_table <- S4Vectors::metadata(deSp_se)[["sig_deSp_result"]]
    significant <- S4Vectors::metadata(deSp_se)[["significant"]]
    p_cutoff <- S4Vectors::metadata(deSp_se)[["p_cutoff"]]
    FC_cutoff <- S4Vectors::metadata(deSp_se)[["FC_cutoff"]]
    ## lollipop
    if (is.character(sig_table)) {
        lipidPlot <- list(
            in.lipidPlot="Insufficient number of significant lipids for plotting lipid plot.",
            lipidPlot="Insufficient number of significant lipids for plotting lipid plot.")
        table_de_lipid <- "Insufficient number of significant lipids for plotting lipid plot."
    } else {
        ## if sig lipids > 20, only show top 10
        if (nrow(sig_table) > 20) {
            sig_table <- sig_table %>%
                dplyr::mutate(plotGroup=ifelse(log2FC >= 0, '+', '-'))
            sig_table_pos <- sig_table[sig_table$plotGroup=="+", ] %>%
                dplyr::arrange(desc(log2FC)) %>% head(10)
            sig_table_neg <- sig_table[sig_table$plotGroup=="-", ] %>%
                dplyr::arrange(log2FC) %>% head(10)
            if (isTRUE("sig_FC" %in% colnames(sig_table)) ){
                table_de_lipid <- rbind(sig_table_pos, sig_table_neg) %>%
                    dplyr::select(-c("plotGroup", "mean_ctrl", "mean_exp"))
            } else {
                table_de_lipid <- rbind(sig_table_pos, sig_table_neg) %>%
                    dplyr::select(-c("plotGroup", "mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp", "method", "statistic"))
            }

        } else {
            if (isTRUE("sig_FC" %in% colnames(sig_table)) ){
                table_de_lipid <- sig_table %>%
                    dplyr::select(-c("mean_ctrl", "mean_exp"))
            } else {
                table_de_lipid <- sig_table %>%
                    dplyr::select(-c("mean_ctrl", "mean_exp", "sd_ctrl", "sd_exp", "method", "statistic"))
            }

        }
        lipidPlot <- .lipidPlot(table_de_lipid, significant)
    }
    ## MA plot & volcano plot
    #if (isTRUE("sig_FC" %in% colnames(sig_table)) ){
    if (sum(group_info$group=="ctrl")==1 && sum(group_info$group=="exp")==1){
        table_maVol <- all_table %>%
            dplyr::select(feature, mean_exp, mean_ctrl, 'M'='log2FC') %>%
            dplyr::mutate(A=(log2(mean_exp) + log2(mean_ctrl))/2) %>%
            dplyr::mutate(negLogP=NA, sig_fc.pval_color=ifelse(
                M>log2(FC_cutoff),'up-regulated', ifelse(-M>log2(FC_cutoff), 'down-regulated', 'none')))
        table_maVol$sig_fc.pval_color <- factor(
            table_maVol$sig_fc.pval_color, levels=c("down-regulated", "none", "up-regulated"))
        volcanoPlot <- list(
            in.volcano="Each group must contain at least two samples to plot the volcano plot.",
            volcano="Each group must contain at least two samples to plot the volcano plot.")
    } else {
        table_maVol <- all_table %>%
            dplyr::select(feature, mean_exp, mean_ctrl, pval, padj, negLog10pval, negLog10padj, 'M'='log2FC') %>%
            dplyr::mutate(A=(log2(mean_exp)+log2(mean_ctrl))/2) %>%
            dplyr::mutate(
                negLogP=get(paste0("negLog10", significant)),
                sig_fc.pval_color= ifelse(
                    M>log2(FC_cutoff) & get(paste0("negLog10", significant))>-log10(p_cutoff),'up-regulated',
                    ifelse(-M>log2(FC_cutoff) & get(paste0("negLog10", significant))>-log10(p_cutoff), 'down-regulated','none'))
            )
        table_maVol$sig_fc.pval_color <- factor(table_maVol$sig_fc.pval_color, levels=c("down-regulated", "none", "up-regulated"))
        volcanoPlot <- .volcanoPlot(all_table, table_maVol, significant)
    }
    if (isFALSE(sum(group_info$group=="ctrl")==1 && sum(group_info$group=="exp")==1) ){
        table_maVol %<>% dplyr::select(-c("negLog10pval", "negLog10padj"))
    }
    maPlot <- .maPlot(table_maVol, significant)

    return(list(
        interactive_de_lipid=lipidPlot$in.lipidPlot, interactive_maPlot=maPlot$in.maPlot,
        interactive_volcanoPlot=volcanoPlot$in.volcano,
        static_de_lipid=lipidPlot$lipidPlot, static_maPlot=maPlot$maPlot,
        static_volcanoPlot=volcanoPlot$volcano,
        table_de_lipid=table_de_lipid, table_ma_volcano=table_maVol))
}

.lipidPlot <- function(table_de_lipid, significant){
    ## all sig. lipid is finite
    suppressWarnings({

        lipidPlot <- ggpubr::ggdotchart(
            table_de_lipid, x="feature", y="log2FC", rotate=TRUE, color="white",
            sorting="descending", add="segments", dot.size=2.5,
            legend.title="log2(FC)", xlab=" ", ylab="log2(Fold Change)",
            legend="right", ggtheme=ggpubr::theme_pubr())
        if (isTRUE("sig_FC" %in% colnames(table_de_lipid)) ) {
            lipidPlot <- lipidPlot +
                ggplot2::geom_point(ggplot2::aes(color=log2FC, size=2.5)) +
                ggplot2::guides(size='none') + ggplot2::labs(colour="log2(FC)") +
                ggplot2::scale_colour_gradient2(
                    low="steelblue", mid="white", high="red", midpoint=0)  +
                ggplot2::scale_y_continuous(
                    breaks=c(-5, -2, -1, 0, 1, 2, 5),
                    labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'), limits=c(-6, 6))
        } else if (nrow(table_de_lipid) == 1) {
            lipidPlot <- lipidPlot +
                ggplot2::geom_point(
                    ggplot2::aes(color=ifelse(
                        get(paste0("negLog10", significant)) > 0, "red", "steelblue"), size=2.5)) +
                ggplot2::guides(size="none") +
                ggplot2::theme(legend.position="none")
        } else {
            lipidPlot <- lipidPlot + ggplot2::geom_point(
                ggplot2::aes(color=get(paste0("negLog10", significant)), size=2.5)) +
                ggplot2::guides(size="none") +
                ggplot2::labs(colour=paste0("-log10(", significant, ")")) +
                ggplot2::scale_colour_gradient2(low="steelblue", mid="white", high="red", midpoint=0)
        }
    })
    in.lipidPlot <- plotly::ggplotly(lipidPlot, tooltip="text")
    hover_table <- table_de_lipid %>% dplyr::arrange(log2FC) %>%
        dplyr::select(feature, log2FC, paste0("negLog10", significant) )
    hover_text <- paste(
        "feature :", hover_table$feature, "<br>", "log2(FC) : ",
        round(hover_table$log2FC, 2), "<br>",
        paste0("-log10(", significant, ") :"),
        round(hover_table[, paste0("negLog10", significant)], 2 ))

    in.lipidPlot[["x"]][["data"]][[1]][["text"]] <- hover_text
    return(list(in.lipidPlot=in.lipidPlot, lipidPlot=lipidPlot))
}

.maPlot <- function(table_maVol, significant){
    in.maPlot <- plotly::plot_ly(
        data=table_maVol, x=~as.numeric(A), y=~as.numeric(M),
        type="scatter", mode="markers", color=~ sig_fc.pval_color,
        colors=c("#4169E1", "#DDDDDD", "#FF4500"), showlegend=TRUE,
        marker=list(size=4), hoverinfo="text", text=~ paste(
            "</br>Lipid:", table_maVol$feature, "</br>A value:",
            round(as.numeric(A), 4), "</br>M value:",
            round(as.numeric(M), 4))) %>%
        plotly::layout(
            xaxis=list(
                title="A = (log<sub>2</sub>(exp)+log<sub>2</sub>(ctrl))/2"),
            yaxis=list(title="M = log<sub>2</sub>(exp)-log<sub>2</sub>(ctrl)"),
            title="MA Plot",
            legend=list(
                title=list(text="log2FC Significant"), orientation='h',
                xanchor="center", x=0.5, y=-0.18))
    maPlot <- ggplot2::ggplot(
        table_maVol, ggplot2::aes(
            x=as.numeric(A), y=as.numeric(M), color=sig_fc.pval_color)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(
            values=c("#4169E1", "#DDDDDD", "#FF4500")) +
        ggplot2::labs(
            color='log2FC Significant', y='M = log2(exp)-log2(ctrl)',
            x="A = (log2(exp)+log2(ctrl))/2", title="MA Plot") +
        ggthemes::theme_hc() +
        ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_vline(xintercept=0)
    return(list(in.maPlot=in.maPlot, maPlot=maPlot))
}

.volcanoPlot <- function(all_table, table_maVol, significant){
    in.volcano <-plotly::plot_ly(
        data=table_maVol, x=~as.numeric(M), y=~as.numeric(negLogP),
        type="scatter", mode="markers", color=~sig_fc.pval_color,
        colors=c("#4169E1", "#DDDDDD", "#FF4500"), marker=list(size=4),
        hoverinfo="text", text=~ paste(
            "</br>Lipid:", table_maVol$feature, "</br>A value:",
            round(as.numeric(A), 4), "</br>M value:",
            round(as.numeric(M), 4), "</br>-log10(p-value):",
            round(-log10(pval), 4), "</br>-log10(padj):",
            round(-log10(padj), 4))) %>%
        plotly::layout(
            xaxis=list(
                title="M = log<sub>2</sub>(exp)-log<sub>2</sub>(ctrl)"),
            yaxis=list(title=paste0("-log<sub>10</sub>(", significant, ")")),
            title="Volcano Plot", legend=list(
                title=list(text="Significant lipid"), orientation='h',
                xanchor="center", x=0.5, y=-0.18))

    volcano <- ggplot2::ggplot(
        table_maVol, ggplot2::aes(
            x=as.numeric(M), y=as.numeric(negLogP),
            color=sig_fc.pval_color)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(
            values=c("#4169E1", "#DDDDDD", "#FF4500")) +
        ggplot2::labs(
            color='Significant lipid', y=paste0('-log10(', significant, ')'),
            x="M = log2(exp)-log2(ctrl)", title="Volcano Plot") +
        ggthemes::theme_hc() +
        ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_vline(xintercept=0)
    return(list(in.volcano=in.volcano, volcano=volcano))
}
