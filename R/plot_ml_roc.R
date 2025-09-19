#' @title plot_ml_roc
#' @description This function generates the overall ROC curve for cross-validations
#' (CVs) with varying feature numbers and the ROC curve for the average CVs based
#' on user-selected feature numbers.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @return Return 2 interactive plots, 2 static plots, and 2 data frames.
#' \enumerate{
#' \item interactive_mean_auc & static_mean_auc: ROC curve plots
#' \item initeractive_roc & static_roc: ROC Curve of average CVs plots
#' \item table_mean_auc_plot: ROC data frame of n features
#' \item table_roc: average ROC curve plot of n features
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- plot_ml_roc(ml_se_sub, feature_num=10)

plot_ml_roc <- function(ml_se, feature_num=10){
    .check_ml_outputSE(ml_se)
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    roc_result <- S4Vectors::metadata(ml_se)$roc_result %>% as.data.frame()
    mean_roc_result <- S4Vectors::metadata(ml_se)$mean_roc_result %>% as.data.frame()
    ## plot table
    plotTab <- .plot_roc_pr_tab(roc_result, mean_roc_result, type="ROC", feature_num)
    ## plot mean auc
    in_auc <- plotly::plot_ly(
        plotTab$auc_tab, x=~(1-specificity), y=~sensitivity)%>%
        plotly::add_lines(
            name=~stats::reorder(feature_num_label, as.numeric(feature_num)),
            hoverinfo="text",
            text=~paste("feature_num :", feature_num, "<br>", feature_num_label1)) %>%
        plotly::add_lines(
            x=0:1, y=0:1, line=list(color='black', width=2, dash='dash'), showlegend=FALSE)%>%
        plotly::layout(
            xaxis=list(title="1-specificity"),
            legend=list(title=list(text="Feature number"), y=0.5, x=1.1, font=list(size=9)),
            title=list(size =15, y=0.99, x=0.1, text="ROC curve"))
    static_auc <- ggplot2::ggplot(
        plotTab$auc_tab,
        ggplot2::aes(
            x=(1-specificity),
            y=sensitivity,color=stats::reorder(feature_num_label, as.numeric(feature_num)))) +
        ggplot2::geom_line() +
        ggplot2::geom_abline(intercept=0, slope=1, linetype="dashed") +
        ggplot2::labs(
            color='Feature number', y='sensitivity', x="1-specificity", title="ROC curve") +
        ggthemes::theme_few()
    ## plot roc
    roc_tab <- plotTab$plot_tab %>% dplyr::filter(cv_fold != "mean")
    roc_mean <- plotTab$plot_tab %>% dplyr::filter(cv_fold == "mean")
    auc_data <- roc_result[c("ranking_method", "ml_method", "cv_fold", "feature_num", "ROC_AUC")] %>%
        dplyr::rename(auc="ROC_AUC")
    auc_stat <- .auc_stats(auc_data, feature=NULL)
    in_roc <- plotly::plot_ly(
        roc_tab, x=~(1-specificity), y=~sensitivity, hoverinfo=NULL)%>%
        plotly::add_trace(
            data=roc_tab, color=~cv_fold, colors ="gray",
            hoverinfo="text", mode="lines", type='scatter',
            text=~paste("cv_fold :", cv_fold))%>%
        plotly::add_trace(
            data=roc_mean, name="mean", x=~(1-specificity), y =~sensitivity,
            line=list(color="red"), hoverinfo="text", mode="lines", type='scatter',
            text=~paste(
                "cv_fold :", cv_fold, "<br>AUC = ",
                paste0(as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95), '-',
                       as.character(auc_stat$upper95), ')'))) %>%
        plotly::add_lines(
            x=0:1, y=0:1, line=list(color='black', width=2, dash='dash'), showlegend=FALSE) %>%
        plotly::layout(
            xaxis=list(title="1-specificity"),
            legend=list(title=list(text="cv_fold"), y=0.5, x=1.1, font=list(size=9)),
            title =list(size=15, y=0.99, x=0.1, text=stringr::str_c(
                'ROC curve for ', as.character(feature_num), ' feature model'))) %>%
        plotly::add_annotations(
            x=0.6, y=0.15, xref="x", yref="y", text=stringr::str_c(
                'AUC=', as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95),
                '-', as.character(auc_stat$upper95), ')\n', 'pvalue=',
                as.character(auc_stat$AUC_pvalue)), xanchor='left', showarrow=FALSE)

    color_break <- unique(plotTab$plot_tab$cv_fold)
    color <- ifelse(color_break=='mean','red','gray')
    static_roc <- ggplot2::ggplot(
        plotTab$plot_tab, ggplot2::aes(x=(1-specificity), y=sensitivity, color=cv_fold)) +
        ggplot2::geom_line() +
        ggplot2::annotate(
            "text", x=0.6, y=0.15, label=stringr::str_c(
                'AUC=', as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95), '-',
                as.character(auc_stat$upper95), ')\n', 'pvalue=', as.character(auc_stat$AUC_pvalue))) +
        ggplot2::scale_color_manual(breaks=color_break, values=color) +
        ggplot2::labs(
            color='cv_fold', y='sensitivity', x="1-specificity",
            title=stringr::str_c('ROC curve for ', as.character(feature_num), ' feature model')) +
        ggthemes::theme_few()

    return(list(
        interactive_mean_auc=in_auc, static_mean_auc=static_auc,
        initeractive_roc=in_roc, static_roc=static_roc,
        table_mean_auc_plot=plotTab$auc_tab, table_roc=plotTab$plot_tab))
}
