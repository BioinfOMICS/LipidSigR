#' @title plot_ml_probability
#' @description This function computes and visualizes the average predicted
#' probabilities for each sample in the testing data across all cross-validation
#' (CV) runs. The function plots the distribution of predicted probabilities for
#' two reference labels and a confusion matrix that includes both the sample count
#' and proportion.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @return Return 1 interactive plots, 2 static plots, and 2 data frames.
#' \enumerate{
#' \item interactive_probability_plot & static_probability_plot: the distribution
#' of predicted probabilities in two reference labels.
#' \item static_confusion_matrix: A confusion matrix composed of sample number and proportion.
#' \item table_probability_plot: table for plotting probability plot.
#' \item table_confusion_matrix: table for plotting confusion matrix.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- plot_ml_probability(ml_se_sub, feature_num=10)
#'
plot_ml_probability <- function(ml_se, feature_num){
    .check_ml_outputSE(ml_se)
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    model_result <- S4Vectors::metadata(ml_se)$model_result %>% as.data.frame()
    model_result %<>% dplyr::group_by(ID, feature_num) %>%
        dplyr::mutate(pred_prob=mean(pred_prob, na.rm=TRUE)) %>%
        dplyr::select(-cv_fold, -pred_label) %>% unique() %>%
        dplyr::mutate(pred_label=ifelse(pred_prob > 0.5, 1, 0)) %>%
        dplyr::filter(feature_num==feature_num)
    plot_tab <- model_result %>%
        dplyr::mutate(true_label=as.factor(true_label)) %>% as.data.frame()
    ## plotting
    in_prob <- plot_tab %>% plotly::plot_ly(
        x=~true_label, y=~pred_prob, type='violin', box=list(visible=FALSE),
        points=FALSE, showlegend=FALSE, color=~true_label, colors=c("#132B43", "#56B1F7")) %>%
        plotly::add_markers(
            x=~jitter(as.numeric(paste(true_label))), y=~pred_prob,
            text=~paste("Actual group :", true_label, "<br>Predicted probability :",
                        round(pred_prob, 2)), hoverinfo="text") %>%
        plotly::layout(
            xaxis=list(zeroline=FALSE, title="Actual group"),
            yaxis=list(zeroline=FALSE, title="Predicted probabilities"),
            title=list(size =15, y=0.99, x=0.1, text="Average sample probability in all CVs"))
    static_prob <- plot_tab %>% ggplot2::ggplot(
        ggplot2::aes(x=true_label, y=pred_prob, color=true_label)) +
        ggplot2::geom_violin(trim=FALSE) +
        ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2)) +
        ggplot2::scale_color_manual(values=c("#132B43", "#56B1F7")) +
        ggthemes::theme_hc() +
        ggplot2::labs(
            y="Predicted probabilities", x="Actual group",
            title="Average sample probability in all CVs") +
        ggplot2::guides(color="none")

    cm <- caret::confusionMatrix(as.factor(model_result$true_label), as.factor(model_result$pred_label))
    cm_tab <- as.data.frame(cm$table) %>% dplyr::group_by(Reference) %>%
        dplyr::mutate(pct=round(Freq/sum(Freq), 2)) %>%
        dplyr::mutate(label=stringr::str_c(as.character(Freq), ' (', as.character(pct), ')'))
    ## plotting
    static_cm <- cm_tab %>% ggplot2::ggplot(
        ggplot2::aes(x=Reference , y=Prediction, fill=pct))+
        ggplot2::geom_tile() + #change 寫法
        ggplot2::scale_fill_gradient(low="white", high="#08306b") +
        ggplot2::theme_bw() +
        ggplot2::labs(x='Actual group', y='Predicted group', title='Confusion matrix')+
        ggplot2::geom_text(ggplot2::aes(label=label,color=pct > 0.5)) +
        ggplot2::scale_color_manual(guide="none", values=c("black", "white")) +
        ggplot2::guides(fill="none") + ggplot2::coord_equal()
    in_cm <- plotly::ggplotly(static_cm)

    #cm_data$pred_prob <- round(cm_data$pred_prob, 3)
    return(list(
        interactive_probability_plot=in_prob, static_probability_plot=static_prob,
        #interactive_confusion_matrix=in_cm,
        static_confusion_matrix=static_cm,
        table_probability_plot=plot_tab, table_confusion_matrix=cm_tab))
}
