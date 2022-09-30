#' @title cor_network
#' @description A correlation network provides users to interrogate the interaction of features in a machine learning model. Users can choose an appropriate feature number according to previous cross-validation results and the features in the best model (based on ROC-AUC+PR-AUC) will be picked up to compute the correlation coefficients between each other.
#' @description The function computes the correlation coefficients according to the results of \code{\link{model_for_net}} and visualizes the correlation network.
#' @param exp_transform_table The output data frame(list[[1]]) of \code{\link{ML_data_process}}.
#' @param lipid_char_table A data frame of lipid characteristics such as name(feature) of lipid, class of lipid, the total length of lipid, and Fatty acid (FA_) related characteristics. NAs are allowed. The name of the first column must be "feature" (lipid species).
#' @param sig_feature The output of \code{\link{model_for_net}} column \bold{'feature'}.
#' @param node_col The output of \code{\link{model_for_net}} column \bold{'importance'}.
#' @param cor_method A character string indicating which correlation coefficient is to be computed. One of "Pearson", "Kendall", or "spearman".
#' @param edge_cutoff A numeric value between 0 and 1. Only the correlation coefficient larger than it will be shown as a line in the plot.
#' @return Return 1 plot.
#' \enumerate{
#' \item visNetwork_plot: the network of feature importance.
#' }
#' @export
#' @examples
#' data("ML_exp_data")
#' data("ML_lipid_char_table")
#' data("ML_condition_table")
#' exp_data <- ML_exp_data
#' lipid_char_table <- ML_lipid_char_table
#' condition_table <- ML_condition_table
#' char_var <- colnames(lipid_char_table)[-1]
#' ML_data <- ML_data_process(exp_data, group_info = condition_table,
#'                            lipid_char_table, char_var[1],
#'                            exclude_var_missing=TRUE, missing_pct_limit=50,
#'                            replace_zero=TRUE, zero2what='min', xmin=0.5,
#'                            replace_NA=TRUE, NA2what='min', ymin=0.5,
#'                            pct_transform=TRUE, data_transform=TRUE,
#'                            trans_type='log', centering=FALSE, scaling=FALSE)
#' ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
#'                       ML_method='Random_forest', split_prop=0.3, nfold=10)
#' model_net <- model_for_net(ML_data[[2]], ML_method='Random_forest',
#'                            varimp_method='Algorithm-based', ML_output[[8]],
#'                            ML_output[[9]], feature_num=10, nsim=5)
#' cor_network(ML_data[[1]], lipid_char_table, model_net[[2]], model_net[[3]],
#'             cor_method='pearson', edge_cutoff=0)
cor_network <- function(exp_transform_table, lipid_char_table,
                        sig_feature, node_col,
                        cor_method, edge_cutoff){
  
  exp_transform_table <- as.data.frame(exp_transform_table)
  lipid_char_table <- as.data.frame(lipid_char_table)
  if(!is(exp_transform_table[,1], 'character')){
    stop("exp_transform_table first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_transform_table)==2){
    if(!is(exp_transform_table[,1], 'character') | sum(class(exp_transform_table[,-1])%in%c("numeric","integer"))!=1){
      stop("exp_transform_table first column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(!is(exp_transform_table[,1], 'character') | sum(vapply(exp_transform_table[,-1], class,character(1))%in%c("numeric","integer"))!=ncol(exp_transform_table[,-1])){
      stop("exp_transform_table first column type must be 'character',others must be 'numeric'")
    }
  }
  if(nrow(exp_transform_table)!=length(unique(exp_transform_table[,1]))){
    stop("exp_transform_table lipids name (features) must be unique")
  }
  if(ncol(exp_transform_table)<61){
    stop("exp_transform_table at least 60 samples.")
  }
  if(nrow(exp_transform_table)<10){
    stop("exp_transform_table number of lipids names (features) must be more than 10.")
  }
  if(sum(!is.na(exp_transform_table[,-1]))==0 | sum(!is.null(exp_transform_table[,-1]))==0){
    stop("exp_transform_table variables can not be all NULL/NA")
  }

  if(!is(lipid_char_table[,1], 'character')){
    stop("lipid_char_table first column must contain a list of lipids names (features).")
  }
  if(tibble::is.tibble(lipid_char_table)){
    if(nrow(lipid_char_table)!=nrow(unique(lipid_char_table[,1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(lipid_char_table)!=length(unique(lipid_char_table[,1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if("class" %in%colnames(lipid_char_table)){
    if(!is(lipid_char_table[,'class'], 'character')){
      stop("lipid_char_table content of column 'class' must be characters")
    }
  }
  if("totallength" %in%colnames(lipid_char_table)){
    if(!class(lipid_char_table[,'totallength'])%in%c("integer","numeric")){
      stop("lipid_char_table content of column 'totallength' must be numeric")
    }
  }
  if("totaldb" %in%colnames(lipid_char_table)){
    if(!class(lipid_char_table[,'totaldb'])%in%c("integer","numeric")){
      stop("lipid_char_table content of column 'totaldb' must be numeric")
    }
  }
  if("totaloh" %in%colnames(lipid_char_table)){
    if(!class(lipid_char_table[,'totaloh'])%in%c("integer","numeric")){
      stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
    }
  }

  rownames(exp_transform_table) <- exp_transform_table$feature
  edge_table <- exp_transform_table %>% dplyr::filter(feature %in% sig_feature) %>%
    dplyr::select(-1) %>% t() %>% stats::cor(method = cor_method, use = 'pairwise.complete.obs') %>%
    as.data.frame() %>% dplyr::mutate(from=rownames(.)) %>%
    tidyr::gather(-from, key='to', value='cor_coef') %>%
    dplyr::filter(abs(cor_coef)>=edge_cutoff)

  node_table <- data.frame(feature=sig_feature, node_col=node_col) %>%
    dplyr::filter(feature%in%c(edge_table$from, edge_table$to))
  if(!is.null(lipid_char_table)){
    node_table <- node_table %>% dplyr::left_join(lipid_char_table, by='feature')
  }
  node_table <- node_table %>%
    dplyr::select(feature, node_col)
  colnames(node_table)[which(colnames(node_table)=="node_col")] <- "color"
  colnames(node_table)[which(colnames(node_table)=="feature")] <- "id"

  node_table <- node_table %>%
    dplyr::mutate(lebel = id,
           shpae = "diamond",
           title = paste0("<p><b>",id))
  node_table <- dplyr::arrange(node_table,color)
  nodes <- node_table
  n_color <- length(unique(node_table$color))
  unique_color <- unique(node_table$color)
  if(min(node_table$color)>0){
    for (i in seq_len(nrow(node_table))) {
      for (a in seq_len(n_color)) {
        if (node_table$color[i]==unique_color[a]){
          node_table$color[i] <- grDevices::colorRampPalette(c("white" , "red"))(n_color)[a]
        }
      }
    }
  }else if(max(node_table$color)<0){
    for (i in seq_len(nrow(node_table))) {
      for (a in seq_len(n_color)) {
        if (node_table$color[i]==unique_color[a]){
          node_table$color[i] <- grDevices::colorRampPalette(c("blue" , "white"))(n_color)[a]
        }
      }
    }
  }else{
    for (i in seq_len(nrow(node_table))) {
      for (a in seq_len(n_color)) {
        if (node_table$color[i]==unique_color[a]){
          node_table$color[i] <- grDevices::colorRampPalette(c("blue","white" , "red"))(n_color)[a]
        }
      }
    }
  }

  edge_table <- edge_table %>%
    dplyr::mutate(width = abs(cor_coef)) %>%
    dplyr::select(from, to, width, cor_coef) %>%
    dplyr::mutate(width = scales::rescale(width, to = c(1, 5)),
           dashes = FALSE,
           smooth = FALSE,
           shadow = FALSE,
           color = ifelse(cor_coef > 0, '#FFDDAA', '#CCBBFF'),
           title = paste0("<p><b>",edge_table$from," vs ",edge_table$to,"</b><br>cor_coef = ",  round(edge_table$cor_coef,5),"</p>"))
  edge_table <-  edge_table[-which(edge_table$from==edge_table$to),]

  color_ledges_n <- c(1,
                      which(nodes$color==unique(nodes$color)[round(length(unique(nodes$color))/4)])[1],
                      which(nodes$color==unique(nodes$color)[round(2*length(unique(nodes$color))/4)])[1],
                      which(nodes$color==unique(nodes$color)[round(3*length(unique(nodes$color))/4)])[1],
                      which(nodes$color==unique(nodes$color)[length(unique(nodes$color))])[1])

  ledges <- data.frame(color = c("white",node_table$color[color_ledges_n]),
                       label = c("Feature\nimportance",
                                 as.character(round(nodes$color[color_ledges_n],3))),
                       shape =c("box",rep("dot",5)),
                       size = c(30,rep(20,5)),
                       font.size =c(35,rep(30,5)))


  visNet <- visNetwork::visNetwork(node_table, edge_table) %>%
    visNetwork::visLayout(randomSeed = 12) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE))%>%
    visNetwork::visEdges(color = list(highlight = "#C62F4B")) %>%
    visNetwork::visPhysics(solver = 'barnesHut',
               stabilization = TRUE,
               barnesHut = list(gravitationalConstant = -1000)) %>%
    visNetwork::visIgraphLayout(layout = "layout_nicely") %>%
    visNetwork::visInteraction(navigationButtons = TRUE)  %>%
    visNetwork::visEvents(dragEnd ="function () {this.setOptions( { physics: false } );}") %>%
    visNetwork::visLegend(useGroups = FALSE,addNodes=ledges,width = 0.15)
  return(visNetwork_plot = visNet)
}
