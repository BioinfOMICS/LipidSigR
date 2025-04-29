#' @title ml_shap
#' @description This function utilizes the Shapley Additive exPlanations (SHAP)
#' method to rank feature importance based on a user-defined number of features.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @param nsim Numeric. A positive integer indicating the times of simulation.
#' Default is \code{5}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("ml_se_sub")
#' shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
ml_shap <- function(ml_se, feature_num=10, nsim=5){
    .check_ml_outputSE(ml_se)
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    if (!is.numeric(nsim) | isFALSE(.check_numeric_range(nsim, 0, NULL)) ) {
        stop("nsim must be a positive integer.")
    }
    best_model <- S4Vectors::metadata(ml_se)$best_model
    best_model_feature <- S4Vectors::metadata(ml_se)$best_model_feature
    transform <- S4Vectors::metadata(ml_se)$transform
    char <- S4Vectors::metadata(ml_se)$char
    ml_method <- S4Vectors::metadata(ml_se)$ml_method

    ml_data <- .ml_process(ml_se, char, transform, type="transpose")
    std1 <- function(x){
        return ((x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
    }
    model <- best_model[[as.character(feature_num)]]
    best_model_feature <- lapply(best_model_feature, function(x) {
        .lipid_name_replace(x, type="replace")
    })
    best_model_feature <- c(
        TRUE, colnames(ml_data)[-1] %in% best_model_feature[[as.character(feature_num)]])
    data <- ml_data[best_model_feature]

    if(ml_method == 'xgboost'){
        shap_values <- SHAPforxgboost::shap.values(
            xgb_model=model,X_train=as.matrix(data[-1]))
        shap_score <- shap_values$shap_score
        shap_long <- SHAPforxgboost::shap.prep(
            xgb_model=model, X_train=as.matrix(data[-1]), top_n=NULL)
    }else if(ml_method == 'Random_forest'){
        pred <- function(object, newdata) {
            stats::predict(object, data=newdata)$predictions
        }
        shap_score <- fastshap::explain(
            model, X=data[-1], nsim=nsim, pred_wrapper=pred)
    }else if(ml_method == 'SVM'){
        pred <- function(object, newdata) {
            attributes(
                stats::predict(object, newdata=newdata, probability=TRUE))$probabilities[,2]
        }
        shap_score <- fastshap::explain(
            model, X=data[-1], nsim=nsim, pred_wrapper=pred)
    }else if(ml_method %in% c('Lasso', 'Ridge', 'ElasticNet')){
        pred <- function(object, newdata) {
            stats::predict(object, newx=as.matrix(newdata), type='response')[, 1]
        }
        shap_score <- fastshap::explain(model, X=data[-1], nsim=nsim, pred_wrapper=pred)
    }
    ## change
    colnames(data) <- .lipid_name_replace(colnames(data), type="revert")
    shap_score <- shap_score %>% as.data.frame()
    colnames(shap_score) <- .lipid_name_replace(colnames(shap_score), type="revert")

    shap_long <- shap_score %>%
        tidyr::gather(key='variable', value='value') %>%
        dplyr::mutate(
            ID=rep(seq_len(nrow(shap_score)), ncol(shap_score)), .before=variable) %>%
        dplyr::mutate(rfvalue=unlist(data[-1]), .after=value) %>%
        dplyr::group_by(variable) %>% dplyr::mutate(stdfvalue=std1(rfvalue), .after=rfvalue) %>%
        dplyr::mutate(mean_value=mean(abs(value)), .after=stdfvalue) %>%
        dplyr::mutate(variable=as.factor(variable)) %>% as.data.frame() %>%
        dplyr::rename(
            c(shapley_value=value, raw_value=rfvalue, normalized_value=stdfvalue, mean_shapley_value=mean_value))
    # final output
    shap_se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(feature_abundance=as.matrix(data[-1]) ),
        rowData=S4Vectors::DataFrame(data[1]),
        colData=S4Vectors::DataFrame(colnames(data)[-1]),
        metadata=list(
            feature_num=feature_num, nsim=nsim,
            shap_result=shap_long, shap_score=as.data.frame(shap_score)) )
    return(shap_se)
}

#' @title plot_ml_shap
#' @description This function is for visualizing the results of Shapley Additive
#' exPlanations (SHAP) results.
#' @param shap_se A SummarizedExperiment object with results computed by \code{\link{ml_shap}}.
#' @return Return 2 interactive plots, 2 static plots, and 2 tables.
#' \enumerate{
#' \item interactive_feature_importance & static_feature_importance: SHAP feature importance plot.
#' \item interactive_summary_plot & static_summary_plot: SHAP value plot.
#' \item table_feature_importance: table for plotting feature importance plot.
#' \item table_summary_plot: table for plotting SHAP value plot.
#' }
#' @export
#' @examples
#' library(SHAPforxgboost)
#' data("ml_se_sub")
#' shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
#' res <- plot_ml_shap(shap_se)
plot_ml_shap <- function(shap_se){
    .check_inputSE(shap_se, metadata_list=c("feature_num", "nsim", "shap_result", "shap_score"))
    abundance <- .extract_df(shap_se, type="abundance") %>%
        tibble::column_to_rownames(var="feature")
    #colnames(abundance) <- .lipid_name_replace(colnames(abundance), type="replace")
    shap_score <- S4Vectors::metadata(shap_se)$shap_score
    #colnames(shap_score) <- .lipid_name_replace(colnames(shap_score), type="replace")
    shap_long <- S4Vectors::metadata(shap_se)$shap_result
    feature_num <- S4Vectors::metadata(shap_se)$feature_num
    nsim <- S4Vectors::metadata(shap_se)$nsim

    topN <- ifelse((feature_num > 10), 10, feature_num)
    ## plotting
    plot_tab <- shap_long[, c("variable", "mean_shapley_value")] %>% as.data.frame() %>%
        unique() %>% dplyr::arrange(dplyr::desc(mean_shapley_value)) %>% .[seq_len(topN),]
    in_mean_shap <- plot_tab %>% plotly::plot_ly(
        x=~mean_shapley_value, y=~stats::reorder(variable, mean_shapley_value), type="bar",
        orientation="h", hoverinfo="text", marker=list(color=~-mean_shapley_value, colorscale="Blues"),
        text=~paste("Feature :", variable, "<br>Mean shapley value :", round(mean_shapley_value, 2))) %>%
        plotly::layout(
            title="SHAP feature importance",
            xaxis=list(title="mean(|Shapley value|)", nticks=15, showline=TRUE, mirror='all'),
            yaxis=list(title=" "))
    static_mean_shap <- plot_tab %>% ggplot2::ggplot(
        ggplot2::aes(y=mean_shapley_value,
                     x=stats::reorder(variable, mean_shapley_value), fill=-mean_shapley_value)) +
        ggplot2::geom_bar(stat="identity") + ggplot2::coord_flip() +
        ggplot2::ggtitle("SHAP feature importance") +
        ggplot2::ylab("mean(|Shapley value|)") + ggplot2::xlab(" ") +
        ggplot2::guides(fill="none") + ggplot2::theme_bw()
    static_all_shap <- SHAPforxgboost::shap.plot.summary.wrap2(
        shap_score=as.data.frame(shap_score), X=as.matrix(abundance), top_n=topN)
    in_all_shap <- plotly::ggplotly(static_all_shap)
    return(list(
        interactive_feature_importance=in_mean_shap,
        static_feature_importance=static_mean_shap,
        interactive_summary_plot=in_all_shap, static_summary_plot=static_all_shap,
        table_feature_importance=plot_tab, table_summary_plot=shap_score))
}

#' @title plot_shap_sample
#' @description This function plots the SHAP feature importance results of user-selected sample.
#' @param shap_se A SummarizedExperiment object with results computed by \code{\link{ml_shap}}.
#' @param sample_id Numeric. The number of samples to display for each feature. Default is \code{10}.
#' @return Return 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item interactive_sample_feature_importance & static_sample_feature_importance:
#' SHAP feature importance plot.
#' \item table_sample_feature_importance: table for plotting SHAP feature importance.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
#' sample_id_list <- unique(S4Vectors::metadata(shap_se)$shap_result$ID)
#' res <- plot_shap_sample(shap_se, sample_id=sample_id_list[10])
plot_shap_sample <- function(shap_se, sample_id=10){
    .check_inputSE(shap_se, metadata_list=c("feature_num", "nsim", "shap_result", "shap_score"))
    shap_long <- S4Vectors::metadata(shap_se)$shap_result
    if(!is.numeric(sample_id) | !sample_id %in% unique(shap_long$ID)){
        stop("sample_id must be one of the ID number in 'ID' column of shap_result.")
    }
    feature_num <- length(unique(shap_long$variable))
    topN <- ifelse((feature_num > 10), 10,  feature_num)
    plot_tab <- shap_long %>% dplyr::filter(ID == sample_id) %>%
        dplyr::mutate(abs_value=abs(shapley_value)) %>%
        dplyr::arrange(dplyr::desc(abs_value)) %>% .[seq_len(topN),]
    in_sample <- plot_tab %>% plotly::plot_ly(
        x=~shapley_value, y=~reorder(variable, shapley_value), orientation="h",
        type="bar", hoverinfo="text",
        marker=list(color=~-shapley_value, colorscale="Blues"),
        text=~paste("Variable :", variable, "<br>Shapley_value :", round(shapley_value, 2))) %>%
        plotly::layout(
            title=stringr::str_c('Feature importance for Sample ', as.character(sample_id)),
            xaxis=list(title="Shapley value", nticks=40, showline=TRUE, mirror='all'),
            yaxis=list(title="", tickfont=list(size=8), zeroline=FALSE, zerolinewidth=0))
    static_sample <- plot_tab %>% ggplot2::ggplot(
        ggplot2::aes(x=shapley_value, y=stats::reorder(variable, shapley_value),
                     fill=-shapley_value)) +
        ggplot2::geom_col() + ggthemes::theme_hc() +
        ggplot2::labs(
            y='', x="Shapley value",
            title=stringr::str_c('Feature importance for Sample ', as.character(sample_id))) +
        ggplot2::guides(fill="none")
    return(list(
        interactive_sample_feature_importance=in_sample,
        static_sample_feature_importance=static_sample,
        table_sample_feature_importance=plot_tab))
}

#' @title plot_shap_force
#' @description The function stacks the SHAP values for each observation and shows
#' how the final output was obtained as a sum of each predictor’s attributions through
#' the force plot. (Randomly plotting a certain portion of the data is optional in
#' case the dataset is large.)
#' @param shap_se A SummarizedExperiment object with results computed by \code{\link{ml_shap}}.
#' @param top_feature Integer. Top number of features to be shown. If top_feature
#' is set to greater than 10, only the top 10 samples will be displayed. Default is \code{10}.
#' @param cluster_method Character. The clustering method. Allowed methods include
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' Default is \code{'ward.D'}.
#' @param group_num Integer. The numbers of groups to be shown in the plots.
#' Default is \code{10}.
#' @return Return 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item interactive_forcePlot & static_forcePlot: SHAP force plot.
#' \item table_forcePlot: table for plotting force plot.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
#' res <- plot_shap_force(shap_se, top_feature=10, cluster_method="ward.D", group_num=10)
plot_shap_force <- function(
        shap_se, top_feature=10, cluster_method="ward.D", group_num=10){
    .check_inputSE(shap_se, metadata_list=c("feature_num", "nsim", "shap_result", "shap_score"))
    shap_score <- S4Vectors::metadata(shap_se)$shap_score
    if(is.null(cluster_method) | isFALSE(cluster_method %in% c(
        "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) ){
        stop('cluster_method must be one of the strings "ward.D", "ward.D2",
         "single", "complete", "average", "mcquitty", "median", "centroid".')
    }
    if (!is.numeric(group_num) | isFALSE(.check_numeric_range(group_num, 0, nrow(shap_score))) ) {
        stop("group_num must be a numeric value between 0 and ", nrow(shap_score)," 1.")
    }
    if (!is.numeric(top_feature) | isFALSE(.check_numeric_range(top_feature, 0, NULL)) ) {
        stop("top_feature must be a positive value.")
    } else {
        if ( top_feature > 10){
            warning("top_feature is set to greater than 10, only the top 10 samples will be displayed")
        }
    }
    topN <- ifelse((top_feature > 10), 10,  top_feature)
    plot_data <- SHAPforxgboost::shap.prep.stack.data(
        shap_contrib=shap_score, top_n=topN, cluster_method=cluster_method,
        n_groups=group_num)
    f_plot <- .forcePlot(plot_data, x_rank='sorted_id')
    plot_data <- plot_data %>% as.data.frame() %>%
        dplyr::select(ID, group, sorted_id, dplyr::everything())
    colnames(plot_data)[seq_len(3)] <- c('Sample ID', 'Group', 'Sorted ID')
    return(list(
        interactive_forcePlot=f_plot$interactive_forcePlot,
        static_forcePlot=f_plot$static_forcePlot,
        table_forcePlot=plot_data))
}

.forcePlot <- function (
        shapobs, id="sorted_id", y_parent_limit=NULL, y_zoomin_limit=NULL,
        x_rank='sorted_id'){
    shapobs_long <- data.table::melt.data.table(
        shapobs, measure.vars=colnames(shapobs)[!colnames(shapobs) %in% c(id, "group", "ID")])
    shapobs_long <- as.data.frame(shapobs_long)

    plotly_colors <- if (dim(shapobs)[2]-1 <= 12) {
        RColorBrewer::brewer.pal(dim(shapobs)[2]-1, "Paired")
    } else {NULL}
    inPlot <- plotly::plot_ly(
        shapobs_long, x=shapobs_long[, which(colnames(shapobs_long)==x_rank)],
        y=~value, type='bar', color=~variable, hoverinfo="text",
        colors=plotly_colors,
        text=~paste(
            "Sorted ID :", shapobs_long[, which(colnames(shapobs_long)==x_rank)],
            "\nShapley values by feature :", round(value, 3), "\nFeature :", variable)) %>%
        plotly::layout(
            barmode="relative", title='SHAP force plot', xaxis=list(title="Sorted ID"),
            yaxis=list(title="Shapley values by feature:\n (Contribution to the base value)"),
            legend=list(title=list(text="Feature"), y=0.5, x=1.1, font=list(size=9)))
    static_plot <- ggplot2::ggplot(
        shapobs_long, ggplot2::aes_string(x=id, y="value", fill="variable")) +
        ggplot2::geom_col(width=1, alpha=0.9) +
        ggplot2::labs(
            fill='Feature', x='Sorted ID',
            y='SHAP values by feature:\n (Contribution to the base value)',
            title='SHAP force plot') +
        ggplot2::geom_hline(yintercept=0, col="gray40") + ggplot2::theme_bw() +
        ggplot2::coord_cartesian(ylim=y_parent_limit)
    if (dim(shapobs)[2]-1<=12){
        static_plot <- static_plot +
            ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(
                dim(shapobs)[2]-1, 'Paired'))
    } else {
        static_plot <- static_plot + ggplot2::scale_fill_viridis_d(option = "D")
        # viridis color scale
    }
    if (!is.null(y_parent_limit) & !is.null(y_zoomin_limit)) {
        warning("Just notice that when parent limit is set, the zoom in axis
              limit won't work, it seems to be a problem of ggforce.\n")
    }
    return(list(
        interactive_forcePlot=inPlot, static_forcePlot=static_plot,
        table_forcePlot=shapobs_long))
}

#' @title plot_shap_dependence
#' @description This function visualizes the SHAP values against the feature
#' values for each variable.
#' @param shap_se A SummarizedExperiment object with results computed by \code{\link{ml_shap}}.
#' @param feature Character. A character string of the feature name the user selects
#' from the "variable" column in \code{shap_result}. (Use \code{\link{extract_summarized_experiment}}
#' to view \code{shap_result}.)
#' @param shap_feature Character. A character string of the feature name the user selects from the
#' "variable" column of \code{shap_result}. (Use \code{\link{extract_summarized_experiment}}
#' to view \code{shap_result}.)
#' @param interaction_index Character. A character string of the feature name the user selects
#' from the "variable" column of \code{shap_result}. (Use \code{\link{extract_summarized_experiment}}
#' to view \code{shap_result}.)
#' @return Return 1 interactive plot, 1 static plot, 1 table.
#' \enumerate{
#' \item interactive_dependence_plot ＆ static_dependence_plot: SHAP dependence plot.
#' \item table_dependence_plot: table for plotting SHAP dependence plot
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' shap_se <- ml_shap(ml_se_sub, feature_num=10, nsim=5)
#' selected_feature <- as.character(unique(S4Vectors::metadata(shap_se)$shap_result$variable))
#' res <- plot_shap_dependence(
#'     shap_se, feature=selected_feature[1], shap_feature=selected_feature[2],
#'     interaction_index=selected_feature[2])
plot_shap_dependence <- function(shap_se, feature, shap_feature, interaction_index){
    .check_inputSE(shap_se, metadata_list=c("feature_num", "nsim", "shap_result", "shap_score"))
    shap_long <- S4Vectors::metadata(shap_se)$shap_result
    if(is.null(feature) | isFALSE(feature %in% unique(shap_long$variable)) ){
        stop("feature must be one of the feature names in the 'variable' column of shap_result.")
    }
    if(is.null(shap_feature) | isFALSE(shap_feature %in% unique(shap_long$variable)) ){
        stop("shap_feature must be one of the feature names in the 'variable' column of shap_result.")
    }
    if(is.null(interaction_index) | isFALSE(interaction_index %in% unique(shap_long$variable)) ){
        stop("interaction_index must be one of the feature names in the 'variable' column of shap_result.")
    }
    plot_tab <- data.frame(
        a=(shap_long %>% dplyr::filter(variable==feature) %>% .$raw_value),
        b=(shap_long %>% dplyr::filter(variable==shap_feature) %>% .$shapley_value),
        c=(shap_long %>% dplyr::filter(variable==interaction_index) %>% .$raw_value))
    in_dependence <- plot_tab %>%
        plotly::plot_ly(x=~a,hoverinfo=NULL) %>%
        plotly::add_trace(
            y=~b, color=~c, hoverinfo="text", mode="markers", type='scatter',
            colors=~colorRampPalette(c("#E6E6FF", "#0000FF"))(length(c)),
            text =~paste(feature, ":", a, "<br>Shapley_value :", round(b, 3),
                         "<br>Feature value :", round(c, 3))) %>%
        plotly::add_lines(
            y=~fitted(loess(b~a)), hoverinfo="none", showlegend=FALSE) %>%
        plotly::colorbar(
            title=stringr::str_c(interaction_index, "\n(Feature value)")) %>%
        plotly::layout(
            title='Shap dependence plot',
            xaxis=list(title=feature, nticks=20, showline=TRUE, mirror='all'),
            yaxis=list(title=stringr::str_c("Shapley value for ", shap_feature),
                       zeroline=FALSE, zerolinewidth=0))
    static_dependence <- plot_tab %>%
        ggplot2::ggplot(ggplot2::aes(x=a, y=b, color=c)) +
        ggplot2::geom_jitter() +
        ggplot2::labs(
            y=stringr::str_c("Shapley value for ", shap_feature), x=feature,
            color=stringr::str_c(interaction_index, "\n(Feature value)")) +
        ggplot2::theme_bw()

    return(list(
        interactive_dependence_plot=in_dependence,
        static_dependence_plot=static_dependence,
        table_dependence_plot=plot_tab))
}

