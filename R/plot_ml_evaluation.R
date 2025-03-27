#' @title plot_ml_evaluation
#' @description This function plots critical indicators to help users evaluate
#' model performance, including the average values and 95% confidence intervals
#' for metrics such as accuracy, sensitivity (recall), specificity, positive
#' predictive value (precision), negative predictive value, F1 score, prevalence,
#' detection rate, detection prevalence, and balanced accuracy.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param eval_method Character. The evaluation method to be used. Allowed methods
#' include 'Accuracy', 'Sensitivity', 'Specificity', 'Pos Pred Value',
#' 'Neg Pred Value', 'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
#' 'Detection Prevalence', 'Balanced Accuracy'. Default is \code{'Accuracy'}.
#' @return Return 1 interactive plot, 1 static plot, and 2 tables.
#' \enumerate{
#' \item interactive_evaluation_plot & static_evaluation_plot: model performance plot
#' \item table_evaluation: table of model evaluation information.
#' \item table_evaluation_plot: table for plotting evaluation plots.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- plot_ml_evaluation(ml_se_sub, eval_method='Accuracy')
plot_ml_evaluation <- function(ml_se, eval_method){
    .check_ml_outputSE(ml_se)
    eval_method_option <- c(
        'Accuracy', 'Sensitivity', 'Specificity', 'Pos Pred Value',
        'Neg Pred Value', 'Precision', 'Recall', 'F1', 'Prevalence',
        'Detection Rate', 'Detection Prevalence', 'Balanced Accuracy')
    if(is.null(eval_method) | isFALSE(eval_method %in% eval_method_option)){
        stop('eval_method must be one of the strings "Accuracy", "Sensitivity",
         "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision",
         "Recall", "F1", "Prevalence", "Detection Rate",
         "Detection Prevalence", or "Balanced Accuracy".')
    }
    cm <- S4Vectors::metadata(ml_se)$confusion_matrix %>% as.data.frame()
    fil_data <- cm %>% dplyr::filter(index==eval_method)
    # evaluation data
    eval_data <- fil_data %>% dplyr::group_by(feature_num) %>%
        dplyr::mutate(mean_value=mean(value)) %>%
        dplyr::arrange(feature_num, cv_fold) %>%
        dplyr::mutate(across(c(value, mean_value), ~round(., 2)))
    # plot table
    plot_tab <- fil_data %>% dplyr::group_by(feature_num) %>%
        dplyr::summarise(mean_value=mean(value), sd_value=sd(value), .groups="drop") %>%
        dplyr::mutate(
            across(c(mean_value, sd_value), ~round(., 2)),
            highest=dplyr::if_else(mean_value==max(mean_value), 'red', 'black'),
            lower=mean_value - sd_value, upper=mean_value + sd_value)
    # plotting
    in_eval <- plot_tab %>% plotly::plot_ly(
        x=~feature_num, y=~mean_value, line=list(color='black'),
        marker=~list(color=highest, size=10), type='scatter', mode='lines',
        showlegend=FALSE, error_y=~list(array=sd_value, color='#000000'),
        text=~paste(
            "feature_num :", feature_num, "<br>mean :", mean_value, "<br>sd :",
            sd_value), hoverinfo="text") %>%
        plotly::add_annotations(
            text=~mean_value, xref="x", yref="y", showarrow=TRUE, arrowhead=4,
            arrowsize=.3, ax=-20, ay=50) %>%
        plotly::layout(
            xaxis=list(title="Number of features", showgrid=TRUE, nticks=20,
                       showline=TRUE, mirror='all'),
            yaxis=list(title="value (%)", range=c(0, 100), nticks=15, mirror='all'),
            title=list(size =15, y=0.99, x=0.7, text=eval_method))
    static_eval <- plot_tab %>% ggplot2::ggplot(
        ggplot2::aes(feature_num, mean_value, ymin=lower, ymax=upper)) +
        ggplot2::geom_errorbar(position="dodge") +
        ggplot2::geom_point(ggplot2::aes(color=highest),size=3) +
        ggplot2::geom_line() +
        ggrepel::geom_text_repel(
            ggplot2::aes(label=mean_value), size=5, box.padding=ggplot2::unit(0.5, "lines"),
            point.padding=ggplot2::unit(0.5, "lines")) +
        ggplot2::scale_color_manual(
            breaks=c('red', 'black'), values=c('red', 'black')) +
        ggplot2::theme_classic(base_size=10) + ggplot2::xlim(0, 100) +
        ggplot2::ylim(0, 100) +
        ggplot2::labs(y="value (%)", x="Number of features", title=eval_method) +
        ggplot2::guides(color="none")

    return(list(
        interactive_evaluation_plot=in_eval, static_evaluation_plot=static_eval,
        table_evaluation=eval_data, table_evaluation_plot=plot_tab))
}
