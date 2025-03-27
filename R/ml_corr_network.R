#' @title ml_corr_network
#' @description Select feature importance from the result of \code{\link{ml_model}}.
#' Two methods are available to rank feature importance: algorithm-based and SHAP
#' analysis. The average importance of the top features across all CV runs will
#' be displayed by specifying a feature count.
#' @param ml_se A SummarizedExperiment object with results computed by \code{\link{ml_model}}.
#' @param feature_importance Character. The method of feature importance. Allowed methods are
#' 'Algorithm-based' and 'SHAP'. Default is \code{'Algorithm-based'}.
#' @param correlation Character. The method for computing correlation coefficient.
#' Allowed methods includes 'pearson', 'kendall', and 'spearman'. Default is \code{'pearson'}.
#' @param edge_cutoff Numeric. A value between 0 and 1. Only the correlation
#' coefficient larger than it will be shown as a line in the plot.
#' @param feature_num Numeric. The number of features to be shown in the plots.
#' A feature number value selected from \code{feature_option} of the
#' SummarizedExperiment object returned by \code{\link{ml_model}}. Usually be
#' one of 2, 3, 5, 10, 20, 50, 100. Default is \code{10}.
#' @param nsim Integer. The times of simulation. Default is \code{5}.
#' @return Return 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item interactive_correlation_network & static_correlation_network: correlation network.
#' \item edge_table & node_table: table for plotting network.
#' }
#' @export
#' @examples
#' data("ml_se_sub")
#' res <- ml_corr_network(ml_se_sub, feature_importance='Algorithm-based',
#'     correlation='pearson', edge_cutoff=0, feature_num=10, nsim=5)
#'
ml_corr_network <- function(
        ml_se, feature_importance='Algorithm-based', correlation='pearson',
        edge_cutoff=0, feature_num=10, nsim=5){
    .check_ml_outputSE(ml_se)
    if(is.null(feature_importance) | isFALSE(feature_importance %in% c("Algorithm-based", "SHAP")) ){
        stop('feature_importance must be "Algorithm-based" or "SHAP".')
    }
    if(is.null(correlation) | isFALSE(correlation %in% c("pearson", "kendall", "spearman")) ){
        stop('correlation must be one of "pearson", "kendall", or "spearman".')
    }
    if (!is.numeric(feature_num) | isFALSE(feature_num %in% S4Vectors::metadata(ml_se)$feature_option) ) {
        stop("feature_num must be a numeric value chosen from the following options: ", paste(S4Vectors::metadata(ml_se)$feature_option, collapse = ", "), " .")
    }
    if (!is.numeric(nsim) | isFALSE(.check_numeric_range(nsim, 0, NULL)) ) {
        stop("nsim must be a integer >= 0.")
    }
    if (!is.numeric(edge_cutoff) | isFALSE(.check_numeric_range(edge_cutoff, 0, 1)) ) {
        stop("edge_cutoff must be a numeric value between 0 and 1.")
    }
    best_model <- S4Vectors::metadata(ml_se)$best_model
    best_model_feature <- S4Vectors::metadata(ml_se)$best_model_feature
    transform <- S4Vectors::metadata(ml_se)$transform
    char <- S4Vectors::metadata(ml_se)$char
    ml_method <- S4Vectors::metadata(ml_se)$ml_method

    if(feature_importance == 'Algorithm-based'){
        num <- which(as.numeric(names(best_model)) == feature_num)
        if (ml_method == 'Random_forest'){
            varimp <- best_model[[num]]$variable.importance
            sig_feature <- names(varimp)
            node_col <- varimp
        }else if(ml_method == 'SVM'){
            varimp <- stats::coef(best_model[[num]])[-1]
            sig_feature <- names(varimp)
            node_col <- varimp
        }else if (ml_method %in% c('Lasso', 'Ridge', 'ElasticNet')) {
            varimp <- stats::coef(best_model[[num]],s = 'lambda.min') %>%
                as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
            sig_feature <- rownames(varimp)
            node_col <- varimp[[1]]
        }else if(ml_method == 'xgboost'){
            varimp <- xgboost::xgb.importance(best_model[[num]]$feature_names, best_model[[num]])
            sig_feature <- varimp$Feature
            node_col <- varimp$Gain
        }
    }else{
        shap_se <- ml_shap(ml_se, feature_num, nsim)
        varImp_result <- S4Vectors::metadata(shap_se)$shap_result
        # varImp_result <- shap(
        #     data=data, best_model=best_model, best_model_feature=best_model_feature,
        #     ml_method=ml_method, feature_n=feature_num, nsim=nsim)$shap_long
        varImp_dir <- varImp_result %>%
            dplyr::mutate(dir=shapley_value*raw_value) %>%
            dplyr::group_by(variable) %>% dplyr::summarise(dir=sum(dir))
        varImp_result <- varImp_result[c(2, 6)] %>% unique() %>%
            dplyr::left_join(varImp_dir,by='variable') %>%
            dplyr::mutate(
                mean_shapley_value=ifelse(dir > 0, mean_shapley_value, -mean_shapley_value))
        sig_feature <- varImp_result$variable
        node_col <- varImp_result$mean_shapley_value
    }
    sig_feature <- .lipid_name_replace(sig_feature, type="revert")

    ml_data <- .ml_process(ml_se, char, transform, type="normal")
    #lipid_char_table <- .extract_df(ml_se, type="lipid")
    rownames(ml_data) <- ml_data$feature

    net_res <- .ml_net(ml_data, sig_feature, node_col, correlation, edge_cutoff)
    return(net_res)
}

.ml_net <- function(ml_data, sig_feature, node_col, correlation, edge_cutoff){
    edge_table <- ml_data %>%
        dplyr::filter(feature %in% sig_feature) %>%
        dplyr::select(-1) %>% t() %>%
        stats::cor(method=correlation, use='pairwise.complete.obs') %>%
        as.data.frame() %>% dplyr::mutate(from=rownames(.)) %>%
        tidyr::gather(-from, key='to', value='cor_coef') %>%
        dplyr::filter(abs(cor_coef) >= edge_cutoff)
    node_table <- data.frame(feature=sig_feature, node_col=node_col) %>%
        dplyr::filter(feature%in%c(edge_table$from, edge_table$to))
    node_table <- node_table %>% dplyr::select(feature, node_col)
    colnames(node_table)[which(colnames(node_table) == "node_col")] <- "color"
    colnames(node_table)[which(colnames(node_table) == "feature")] <- "id"
    node_table <- node_table %>%
        dplyr::mutate(lebel=id, shpae="diamond", title=paste0("<p><b>", id))
    node_table <- dplyr::arrange(node_table, color)
    nodes <- node_table
    n_color <- length(unique(node_table$color))
    unique_color <- unique(node_table$color)
    node_table$color_label <- node_table$color

    for (i in seq_len(nrow(node_table))) {
        for (a in seq_len(n_color)) {
            if (min(node_table$color) > 0) {
                color_label <- grDevices::colorRampPalette(c("white", "red"))(n_color)[a]
            } else if(max(node_table$color) < 0) {
                color_label <- grDevices::colorRampPalette( c("blue", "white"))(n_color)[a]
            } else {
                color_label <- grDevices::colorRampPalette(c("blue", "white", "red"))(n_color)[a]
            }
            if (node_table$color[i] == unique_color[a]){
                node_table$color[i] <- color_label
            }
        }
    }
    edge_table <- edge_table %>%
        dplyr::mutate(width=abs(cor_coef)) %>%
        dplyr::select(from, to, width, cor_coef) %>%
        dplyr::mutate(
            width=scales::rescale(width, to=c(1, 5)), dashes=FALSE, smooth=FALSE,
            shadow=FALSE, color=ifelse(cor_coef > 0, '#FFDDAA', '#CCBBFF'),
            title=paste0(
                "<p><b>", edge_table$from, " vs ", edge_table$to, "</b><br>cor_coef = ",
                round(edge_table$cor_coef, 5), "</p>"))
    edge_table <-  edge_table[-which(edge_table$from == edge_table$to),]
    net_res <- .ml_net_plot(nodes, node_table, edge_table, unique_color)
    return(net_res)
}

.ml_net_plot <- function(nodes, node_table, edge_table, unique_color){
    ## plotting
    color_ledges_n <- c(
        1, which(nodes$color==unique(nodes$color)[round(length(unique(nodes$color))/4)])[1],
        which(nodes$color==unique(nodes$color)[round(2*length(unique(nodes$color))/4)])[1],
        which(nodes$color==unique(nodes$color)[round(3*length(unique(nodes$color))/4)])[1],
        which(nodes$color==unique(nodes$color)[length(unique(nodes$color))])[1])
    ledges <- data.frame(
        color=c("white", node_table$color[color_ledges_n]),
        label=c("Feature\nimportance", as.character(round(nodes$color[color_ledges_n], 3))),
        shape=c("box", rep("dot", 5)), size=c(30, rep(20, 5)), font.size=c(35, rep(30, 5)))
    in_net <- visNetwork::visNetwork(node_table, edge_table) %>%
        visNetwork::visLayout(randomSeed=12) %>%
        visNetwork::visOptions(
            highlightNearest=list(enabled=TRUE, degree=1, hover=TRUE))%>%
        visNetwork::visEdges(color=list(highlight="#C62F4B")) %>%
        visNetwork::visPhysics(
            solver='barnesHut', stabilization=TRUE, barnesHut=list(gravitationalConstant=-1000)) %>%
        visNetwork::visIgraphLayout(layout="layout_nicely") %>%
        visNetwork::visInteraction(navigationButtons=TRUE)  %>%
        visNetwork::visEvents(
            dragEnd="function(){this.setOptions({physics:false});}") %>%
        visNetwork::visLegend(useGroups=FALSE, addNodes=ledges, width=0.15)

    unique_color_label <- round(unique_color,3)
    unique_color <- unique(node_table$color)
    edge_color <- unique(edge_table$color)
    colnames(node_table)[which(colnames(node_table)=='color')] <- "Feature importance"
    ig <- igraph::graph_from_data_frame(
        d=edge_table, vertices=node_table, directed=FALSE)
    tg <- tidygraph::as_tbl_graph(ig) %>% tidygraph::activate(nodes) %>%
        dplyr::mutate(label=name)
    static_net <- tg %>% ggraph::ggraph(layout="fr") +
        ggraph::geom_edge_arc(
            colour="gray50", lineend="round", strength=.1, alpha=.1) +
        ggraph::geom_edge_fan(ggplot2::aes(color=color), show.legend=FALSE) +
        ggraph::scale_edge_color_manual(breaks=edge_color, values=edge_color) +
        ggraph::geom_node_point(size=6, ggplot2::aes(color=`Feature importance`)) +
        ggplot2::scale_color_manual(
            breaks=unique_color, values=unique_color, labels=unique_color_label) +
        ggraph::geom_node_text(
            ggplot2::aes(label=label), repel=TRUE,
            point.padding=ggplot2::unit(0.2, "lines"), colour="gray10") +
        ggraph::theme_graph(background="white") +
        ggplot2::guides(edge_width="none", edge_alpha="none") +
        ggraph::theme_graph(base_family="sans")
    return(list(
        interactive_correlation_network=in_net,
        static_correlation_network=static_net,
        edge_table=edge_table, node_table=node_table))
}
