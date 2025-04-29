#' @title plot_ml_feature
#' @description This function computes and ranks each feature's contribution
#' based on the user-defined number of features and visualizes their importance.
#' When users select a specific number of features, the frequency of the top
#' features selected across all CV runs is displayed.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @return Return 2 interactive plots, 2 static plots, and 2 data frames.
#' \enumerate{
#' \item interactive_selected_frequency & static_selected_frequency: selected frequency plot
#' \item interactive_feature_importance & static_feature_importance: feature importance plot
#' \item table_selected_frequency: table for plotting selected frequency plot
#' \item table_feature_importance: table for plotting feature importance plot
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- plot_ml_feature(ml_se_sub, feature_num=10)
plot_ml_feature <- function(ml_se, feature_num=10){
    .check_ml_outputSE(ml_se)
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    selected_features <- S4Vectors::metadata(ml_se)$selected_features
    varimp_result <- S4Vectors::metadata(ml_se)$feature_importance_result %>% as.data.frame()
    ranking_method <- S4Vectors::metadata(ml_se)$ranking_method
    ml_method <- S4Vectors::metadata(ml_se)$ml_method
    nfold <- S4Vectors::metadata(ml_se)$nfold

    if (is.null(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }

    feature_loc <- which(names(selected_features[[1]])==as.character(feature_num))
    cv_sele_feature <- selected_features %>%
        purrr::map(function(x){x[[feature_loc]]}) %>% unlist %>% table()
    freq_data <- data.frame(
        feature=names(cv_sele_feature), sele_freq=as.numeric(cv_sele_feature)/nfold) %>%
        dplyr::mutate(feature_num=feature_num, .after=feature) %>%
        dplyr::mutate(ml_method=unique(varimp_result$ml_method), .before=feature_num) %>%
        dplyr::mutate(ranking_method=unique(varimp_result$ranking_method), .before=ml_method) %>%
        dplyr::arrange(dplyr::desc(sele_freq))
    imp_data <- varimp_result %>% dplyr::filter(feature_num==feature_num) %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(importance=mean(importance, na.rm=TRUE)) %>%
        dplyr::select(-cv_fold) %>% dplyr::distinct() %>%
        dplyr::arrange(dplyr::desc(importance))
    if (nrow(freq_data)>=10) {
        plot_freq <- freq_data[seq_len(10), ]
        plot_imp <- imp_data[seq_len(10), ]
    } else {
        plot_freq <- freq_data
        plot_imp <- imp_data
    }
    imp_data$importance <- round(imp_data$importance, 3)
    ## plotting
    in_freq <- plot_freq %>% plotly::plot_ly(
            x=~(sele_freq*100), y=~stats::reorder(feature, sele_freq), type="bar",
            orientation="h", marker=list(color=~-sele_freq, colorscale="Blues"),
            text=~paste("Feature :", feature, "<br>Selected frequency :", sele_freq*100, "%"),
            hoverinfo="text") %>%
        plotly::layout(
            xaxis=list(title="Selected frequency (%)", showgrid=TRUE, nticks=20,
                showline=TRUE, mirror='all'),
            yaxis=list(title=""),
            title =list(size=15, y=0.99, x=0.1, text=stringr::str_c(
                'Selected frequency for ', as.character(feature_num), ' feature model')) )
    static_freq <- plot_freq %>% ggplot2::ggplot(
        ggplot2::aes(x=(sele_freq*100), y=stats::reorder(feature, sele_freq),
                     fill=-sele_freq)) +
        ggplot2::geom_col() + ggplot2::labs(y='', x="Selected frequency (%)") +
        ggplot2::guides(fill="none") + ggthemes::theme_hc()
    ## plotting
    in_imp <- plot_imp %>% plotly::plot_ly(
        x=~importance,y=~stats::reorder(feature, importance), type="bar",
        orientation="h", marker=list(color=~-importance, colorscale="Blues"),
        text=~paste("Feature :", feature, "<br>Average importance :",
                    round(importance, 2)), hoverinfo="text") %>%
        plotly::layout(
            xaxis=list(title='Average importance', showgrid=TRUE, nticks=20,
                       showline=TRUE, mirror='all'),
            yaxis=list(title=""),
            title =list(size=15, y=0.99, x=0.1, text=stringr::str_c(
                'Feature importance for ', as.character(feature_num), ' feature model')) )
    static_imp <- plot_imp %>% ggplot2::ggplot(
        ggplot2::aes(x=importance, y=stats::reorder(feature, importance), fill=-importance)) +
        ggplot2::geom_col() + ggplot2::labs(y='', x="Average importance") +
        ggplot2::guides(fill="none") + ggthemes::theme_hc()
    return(list(
        interactive_selected_frequency=in_freq, static_selected_frequency=static_freq,
        interactive_feature_importance=in_imp, static_feature_importance=static_imp,
        table_selected_frequency=freq_data, table_feature_importance=imp_data))
}
