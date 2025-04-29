#' @title plot_ml_pr
#' @description This function generates the overall Precision-Recall (PR) curve
#' for cross-validations (CVs) with different feature numbers and the PR curve
#' for the average CVs based on user-selected feature numbers.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @return Return 2 interactive plots, 2 static plots, and 2 data frames.
#' \enumerate{
#' \item interactive_mean_auc & static_mean_auc: PR curve plots.
#' \item initeractive_pr & static_pr: PR Curve of average CVs plots.
#' \item table_mean_auc_plot: data frame of the AUC, recall, and precision of PR of n features.
#' \item table_pr: average PR curve plot of n feature.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- plot_ml_pr(ml_se_sub, feature_num=10)

plot_ml_pr <- function(ml_se, feature_num){
    .check_ml_outputSE(ml_se)
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    pr_result <- S4Vectors::metadata(ml_se)$pr_result %>% as.data.frame()
    mean_roc_result <- S4Vectors::metadata(ml_se)$mean_roc_result %>% as.data.frame()
    plotTab <- .plot_roc_pr_tab(pr_result, mean_roc_result, type="PR", feature_num)
    ## mean auc plot
    in_auc <- plotly::plot_ly(
        plotTab$auc_tab, x=~recall, y=~precision) %>%
        plotly::add_lines(
            name=~stats::reorder(feature_num_label, as.numeric(feature_num)),
            hoverinfo="text", text=~paste(
                "feature_num :", feature_num, "<br>", feature_num_label1)) %>%
        plotly::layout(
            legend=list(title=list(text="Feature number"),
                        y=0.5, x=1.1, font=list(size=9)),
            title=list(size=15, y=0.99, x=0.1, text="PR curve"))
    static_auc <- ggplot2::ggplot(
        plotTab$auc_tab,
        ggplot2::aes(
            x=recall, y=precision, color=stats::reorder(
                feature_num_label,as.numeric(feature_num)))) +
        ggplot2::geom_line() +
        ggplot2::labs(
            color='Feature number', y='precision', x="recall", title="PR curve") +
        ggthemes::theme_few()
    ## pr plot
    pr_tab <- plotTab$plot_tab %>% dplyr::filter(cv_fold != "mean")
    pr_mean <- plotTab$plot_tab %>% dplyr::filter(cv_fold == "mean")
    auc_data <- pr_result[c("ranking_method", "ml_method", "cv_fold", "feature_num", "PR_AUC")] %>%
        dplyr::rename(auc="PR_AUC")
    auc_stat <- .auc_stats(auc_data, feature=NULL)
    in_pr <- plotly::plot_ly(
        pr_tab, x=~recall, y=~precision, hoverinfo=NULL) %>%
        plotly::add_lines(
            color=~cv_fold, colors="gray", hoverinfo="text", text=~paste("cv_fold :", cv_fold))%>%
        plotly::add_lines(
            data=pr_mean, name="mean", x=~recall, y=~precision,
            line=list(color="red"), hoverinfo="text",
            text=~paste(
                "cv_fold :", cv_fold, "<br>",
                paste0(as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95), '-',
                       as.character(auc_stat$upper95), ')'))) %>%
        plotly::layout(
            legend=list(title=list(text="cv_fold"), y=0.5, x=1.1, font=list(size=12)),
            title=list(size=15, y=0.99, x=0.1, text=stringr::str_c(
                'PR curve for ', as.character(feature_num) , ' feature model'))) %>%
        plotly::add_annotations(
            x=0.2, y=0.3, xref="x", yref="y", text=stringr::str_c(
                'AUC=', as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95), '-',
                as.character(auc_stat$upper95), ')\n', 'pvalue=', as.character(auc_stat$AUC_pvalue)),
            xanchor='left', showarrow=FALSE)

    color_break <- unique(plotTab$plot_tab$cv_fold)
    color <- ifelse(color_break=='mean','red','gray')
    static_pr <- ggplot2::ggplot(
        plotTab$plot_tab, ggplot2::aes(x=recall, y=precision, color=cv_fold)) +
        ggplot2::geom_line() +
        ggplot2::annotate(
            "text", x=0.2, y=0.3, label=stringr::str_c(
                'AUC=', as.character(auc_stat$mean_AUC), '(', as.character(auc_stat$lower95), '-',
                as.character(auc_stat$upper95), ')\n', 'pvalue=', as.character(auc_stat$AUC_pvalue))) +
        ggplot2::scale_color_manual(breaks=color_break, values=color) +
        ggplot2::labs(
            color='cv_fold', y='precision', x="recall", title=stringr::str_c(
                'PR curve for ', as.character(feature_num) , ' feature model')) +
        ggthemes::theme_few()

    return(list(
        interactive_mean_auc=in_auc, static_mean_auc=static_auc,
        initeractive_pr=in_pr, static_pr=static_pr,
        table_mean_auc_plot=plotTab$auc_tab, table_pr=plotTab$plot_tab))
}
