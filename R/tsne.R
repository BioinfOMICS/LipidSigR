#' @title tsne
#' @description A method of dimensionality reduction to transform data from a
#' high-dimensional space into a low-dimensional space so that the vital
#' properties of the original data can be retained.
#' @param exp_transform_table A data frame of predictors, including
#' features (molecule, lipid class, etc.) and their expression of each sample.
#' NAs are not allowed. The name of the first column must be
#' "feature" (lipid species).
#' @param group_info A data frame with sample name, label name, the group name
#' of samples (ctrl & exp), and pair number representing 'the pair' for the
#' t-test/Wilcoxon test. NAs are allowed.
#' @param sig_feature A character string for identifying significant variable
#' in \bold{exp_transform_table}.
#' @param pca logical; Whether an initial PCA step should be
#' performed (default: TRUE).
#' @param perplexity numeric; Perplexity parameter
#' (should not be bigger than 3 * perplexity < nrow(X) - 1.
#' @param max_iter A integer; Number of iterations.
#' @param cluster_method A character string indicating which method to be used
#' for clustering. Allowed method include \bold{"kmeans"}, \bold{"kmedoids"},
#' \bold{"hclustering"}, \bold{"dbscan"}, \bold{"group_info"}.
#' @param group_num A positive integer specifying the number of clusters.
#' Group number must be between 1 and 10.
#' @param var1 A character string or numeric.
#' According to the selected cluster method.
#' \itemize{
#'   \item \emph{kmeans}: specifying the metric to be used for calculating
#'   dissimilarities between observations include \bold{"euclidean"}
#'   and \bold{"manhattan"}.
#'   \item \emph{hclustering}: which correlation coefficient is to be computed
#'   or which the distance measure to be used include \bold{"pearson"},
#'   \bold{"kendall"}, \bold{"spearman"}, \bold{"euclidean"},
#'   \bold{"manhattan"}, \bold{"maximum"}, \bold{"canberra"}, \bold{"binary"},
#'   \bold{"minkowski"}.
#'   \item \emph{dbscan}: size of the epsilon neighborhood.
#' }
#' @param var2 A character string or numeric. According to the selected
#' cluster method.
#' \itemize{
#'   \item \emph{hclustering}: The agglomeration method to be used. This
#'   should be (an unambiguous abbreviation of) one of \bold{"ward.D"},
#'   \bold{"ward.D2"}, \bold{"single"}, \bold{"complete"},
#'   \bold{"average"}(=UPGMA), \bold{"mcquitty"}(= WPGMA),
#'   \bold{"median"}(= WPGMC), or \bold{"centroid"} (= UPGMC).
#'   \item \emph{dbscan}: number of minimum points in the
#'   eps region (for core points).
#' }
#' @param insert_ref_group A character string. The name of 'ctrl'
#' after name conversion.
#' @param ref_group A character string. The name of 'exp'
#' after name conversion.
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @return Return a list with 1 data frame and 1 plot.
#' \enumerate{
#' \item tsne_data: A data frame of tsne data.
#' \item tsne_plot: tsne plot.
#' }
#' @export
#' @examples
#' data("profiling_exp_data")
#' exp_data <- profiling_exp_data
#' exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
#'                                     missing_pct_limit=50,
#'                                     replace_zero=TRUE, zero2what='min',
#'                                     xmin=0.5, replace_NA=TRUE,
#'                                     NA2what='min', ymin=0.5,
#'                                     pct_transform=TRUE, data_transform=TRUE,
#'                                     trans_type='log', centering=FALSE,
#'                                     scaling=FALSE)
#' tsne(exp_transform_table, group_info=NULL, sig_feature=NULL, pca=TRUE,
#'      perplexity=5, max_iter=500, cluster_method='kmeans', group_num=2,
#'      var1='euclidean', var2=NULL, insert_ref_group=NULL,
#'      ref_group=NULL, plotly=TRUE)
tsne <- function(exp_transform_table, group_info=NULL, sig_feature=NULL,
                 pca=TRUE, perplexity=5, max_iter=500, cluster_method,
                 group_num=NULL, var1=NULL, var2=NULL,
                 insert_ref_group=NULL, ref_group=NULL, plotly=TRUE){

  exp_transform_table <- as.data.frame(exp_transform_table)
  if(!is.null(group_info)){
    group_info <- as.data.frame(group_info)
  }
  if(ncol(exp_transform_table) < 3){
    stop("At least 2 samples.")
  }
  if(nrow(exp_transform_table) < 5){
    stop("At least 4 features.")
  }
  if(sum(vapply(exp_transform_table[, -1],
                class,character(1)) %in% c("numeric", "integer")) !=
     ncol(exp_transform_table[, -1])){
    stop("First column type must be 'character',others must be 'numeric'")
  }
  if(tibble::is_tibble(exp_transform_table)){
    if(nrow(exp_transform_table) != nrow(unique(exp_transform_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(exp_transform_table) != length(unique(exp_transform_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if(!is.null(sig_feature)){
    if(length(sig_feature) < 2){
      stop("At least 2 significant lipid features.")
    }
  }
  if(sum(is.na(exp_transform_table[, 1])) > 0){
    stop("The lipids name (features) can not contains NA")
  }
  if(!is.null(group_info)){
    if(ncol(group_info) == 4){
      if(sum(vapply(group_info[, seq_len(3)], class,
                    character(1)) != "character") == 0){
        if("pair" %in% colnames(group_info)){
          if(which(colnames(group_info) == "pair") != 4){
            stop("group_info column must arrange in order of sample_name,
                 label_name, group, pair(optional).")
          }
        }else{
          stop("group_info column must arrange in order of sample_name,
               label_name, group, pair(optional).")
        }
      }else{
        stop("group_info first 3 columns must be characters.")
      }
      if(sum(!is.na(group_info[, 4])) != 0 |
         sum(table(group_info[, 4]) != 2) != 0 &
         sum(is.na(group_info[, 4])) != 0){
        stop("group_info each pair must have a specific number, staring from
             1 to N. Cannot have NA, blank, or skip numbers.")
      }
      if(sum(group_info[, 1] %in%
             colnames(exp_transform_table)) != nrow(group_info) |
         sum(group_info[, 1] %in%
             colnames(exp_transform_table)) != ncol(exp_transform_table[, -1])){
        stop("group_info 'sample_name' must same as the name of samples of
             exp_transform_table")
      }
      if(length(unique(group_info[, 3])) == 2){
        if(sum(table(group_info[, 3]) >= 1) != 2){
          stop("group_info column 'group' only can have 2 groups,
               and >= 1 sample for each group.")
        }
      }else{
        stop("group_info column 'group' only can have 2 groups,
             and >= 1 sample for each group.")
      }
      if(!is.null(insert_ref_group)){
        if(!insert_ref_group %in% group_info[, 3]){
          stop("The reference group entered by users must be included in
               the group_info.")
        }
      }
    }else if(ncol(group_info) == 3){
      if("pair" %in% colnames(group_info)){
        stop("group_info column must arrange in order of sample_name,
             label_name, group, pair(optional).")
      }
      if(sum(vapply(group_info,class,character(1)) != "character") != 0){
        stop("group_info first 3 columns must be characters.")
      }
      if(sum(group_info[, 1] %in%
             colnames(exp_transform_table)) != nrow(group_info) |
         sum(group_info[, 1] %in%
             colnames(exp_transform_table)) != ncol(exp_transform_table[,-1])){
        stop("group_info 'sample_name' must same as the name of samples of
             exp_transform_table")
      }
      if(length(unique(group_info[, 3])) == 2){
        if(sum(table(group_info[, 3]) >= 1) != 2){
          stop("group_info column 'group' only can have 2 groups,
               and >= 1 sample for each group.")
        }
      }else{
        stop("group_info column 'group' only can have 2 groups,
             and >= 1 sample for each group.")
      }
      if(!is.null(insert_ref_group)){
        if(!insert_ref_group %in% group_info[, 3]){
          stop("The insert_ref_group entered by users must be included in
               the group_info.")
        }
      }
    }
  }
  if(cluster_method == 'kmeans'){
    if(!is.null(group_num)){
      if( group_num < 1 ||  group_num > 10){
        stop('Group number must be between 1 and 10')
      }
    }else{
      stop('Group number can not be NULL')
    }
  }else if(cluster_method == 'kmedoids'){
    if(!var1 %in% c("euclidean","manhattan")){
      stop('var1 must be one of the strings "euclidean" or "manhattan".')
    }
  }else if(cluster_method == 'hclustering'){
    if(!var1 %in% c('pearson', 'spearman', 'kendall', "euclidean", "maximum",
                    "manhattan", "canberra", "binary", "minkowski")){
      stop('var1 must be one of the strings "pearson", "kendall", "spearman",
           "euclidean", "maximum", "manhattan", "canberra", "binary",
           or "minkowski"')
    }
    if(!var2 %in% c("ward.D", "ward.D2", "single", "complete", "average",
                    "mcquitty", "median", "centroid")){
      stop('var2 must be one of the strings "ward.D","ward.D2","single",
           "complete","average","mcquitty","median","centroid".')
    }
  }else if(cluster_method == 'dbscan'){
    if(!class(var1) %in% c("numeric", "integer")){
      stop("var1 must be 'numeric'")
    }
    if(!class(var2) %in% c("numeric", "integer")){
      stop("var2 must be 'numeric'")
    }
  }

  if(!is.null(insert_ref_group) & !is.null(ref_group) & !is.null(group_info)){
    exp_raw_name <- ref_group[-which(insert_ref_group %in% ref_group)]
    group_info$group[which(group_info$group == 'ctrl')] <-  insert_ref_group
    group_info$group[which(group_info$group == 'exp')] <-  exp_raw_name
  }

  if(!is.null(sig_feature)){
    exp_transform_table <- exp_transform_table %>%
      dplyr::filter(eval(parse(text=colnames(exp_transform_table)[1])) %in%
                      sig_feature)
  }

  num <- apply(exp_transform_table[-1], 1, FUN=function(x){length(unique(x))})
  exp_transform_table <- exp_transform_table[(num!=1),]
  exp_transform_table <- exp_transform_table[!is.infinite(
    rowSums(exp_transform_table[-1], na.rm=TRUE)), ]

  tsne_table <- exp_transform_table[-1] %>% t() %>% as.data.frame()

  colnames(tsne_table) <- exp_transform_table[[1]]


  tsne_table <- stats::na.omit(tsne_table)
  if(nrow(tsne_table) <= (3*perplexity)){
    stop('Samples must be larger than 3*perplexity')
  }


  tsne <- Rtsne::Rtsne(tsne_table, check_duplicates=FALSE , pca=pca,
                       perplexity=perplexity, verbose=TRUE,
                       max_iter=max_iter,theta=0)
  tsne_data <- as.data.frame(tsne$Y) %>%
    dplyr::mutate(sample_name=colnames(exp_transform_table)[-1])
  colnames(tsne_data)[seq_len(2)] <- c('tsne1', 'tsne2')

  #grouping
  if(cluster_method == 'kmeans'){
    cluster_group <- stats::kmeans(tsne$Y, centers=group_num)$cluster
  }else if(cluster_method == 'kmedoids'){
    cluster_group <- cluster::pam(tsne$Y, k=group_num,
                                  metric=var1)$cluster #euclidean manhattan
  }else if(cluster_method == 'hclustering'){
    if(var1 %in% c('pearson', 'spearman', 'kendall')){
      dist.fun <- function(x){
        x <- t(x)
        cor.mat <- stats::cor(x, method=var1, use='complete.obs')
        cor.mat <- (1-cor.mat)
        cor.dist <- stats::as.dist(cor.mat)
        return(cor.dist)
      }
    }else{
      dist.fun <- function(x) stats::dist(x, method=var1)
    }
    cluster_group <- stats::hclust(dist.fun(tsne$Y), method=var2)
    cluster_group <- stats::cutree(cluster_group, k=group_num)
  }else if(cluster_method == 'dbscan'){
    cluster_group <- dbscan::dbscan(tsne$Y, eps=var1, minPts=var2)$cluster
    cluster_group <- ifelse(cluster_group > 9, 9, cluster_group) %>%
      as.character()
    cluster_group <- ifelse(cluster_group == '0', 'noise', cluster_group)
  }else if(cluster_method == 'group_info'){
    group_order <- purrr::map_dbl(rownames(tsne_table),
                                  function(x){which(x ==
                                                      group_info$sample_name)})
    cluster_group <- group_info[group_order, ] %>% .$group
    group_num <- length(unique(cluster_group))
  }

  cluster_group <- as.character(cluster_group)
  tsne_data <- tsne_data %>% dplyr::mutate(group=cluster_group) %>%
    dplyr::select(sample_name, group, dplyr::everything())

  #tsne plot*5
  color <- c("#00AFBB", "#E7B800", "#FC4E07", "#42B540FF", "#BB3099",
             "#EE0099", "#0000AC", "#868686FF", '#00468BFF', 'black')
  if(cluster_method %in% c('dbscan')){
    tsne_plot <- factoextra::fviz_cluster(
      list(data=as.data.frame(tsne$Y), cluster=cluster_group),
      palette=color[seq_len(length(unique(cluster_group)))],
      geom="point",  stand=FALSE, ellipse=FALSE, ggtheme=ggplot2::theme_bw(),
      show.clust.cent=FALSE, xlab='tsne-1', ylab='tsne-2', main="t-SNE",
      legend.title="Groups")
    if(plotly == TRUE){
      tsne_plot <- plotly::ggplotly(tsne_plot)
      for (i in seq_len(length(tsne_plot$x$data))){
        if (!is.null(tsne_plot$x$data[[i]]$name)){
          tsne_plot$x$data[[i]]$name <-  gsub(
            "\\(", "",
            stringr::str_split(tsne_plot$x$data[[i]]$name, ",")[[1]][1])
        }
        if(i<=length(unique(cluster_group))){
          tsne_plot$x$data[[i]]$text <- paste(
            "x :", round(tsne_plot$x$data[[i]]$x, 3),
            "\ny :", round(tsne_plot$x$data[[i]]$y, 3),
            "\nGroups :",tsne_plot$x$data[[i]]$name,
            "\nSample name :", tsne_data[which(
              tsne_data$group == tsne_plot$x$data[[i]]$name), ]$sample_name)
        }
      }
    }

  }else{
    tsne_plot <- factoextra::fviz_cluster(
      list(data=as.data.frame(tsne$Y), sample_id=tsne_data$sample_name,
           cluster=cluster_group),
      palette=color[seq_len(group_num)],
      geom="point", stand=FALSE, ellipse=TRUE, ellipse.type='norm',
      ggtheme=ggplot2::theme_bw(), show.clust.cent=FALSE, xlab='tsne-1',
      ylab='tsne-2', main="t-SNE", legend.title="Groups")

    if(plotly == TRUE){
      tsne_plot <- plotly::ggplotly(tsne_plot)
      for (i in seq_len(length(tsne_plot$x$data))){
        if (!is.null(tsne_plot$x$data[[i]]$name)){
          tsne_plot$x$data[[i]]$name <- gsub(
            "\\(", "",
            stringr::str_split(tsne_plot$x$data[[i]]$name, ",")[[1]][1])
        }
        if(i <= group_num){
          tsne_plot$x$data[[i]]$text <- paste(
            "x :", round(tsne_plot$x$data[[i]]$x, 3),
            "\ny :", round(tsne_plot$x$data[[i]]$y, 3),
            "\nGroups :", tsne_plot$x$data[[i]]$name,
            "\nSample name :", tsne_data[which(
              tsne_data$group == tsne_plot$x$data[[i]]$name), ]$sample_name)
        }else if(i >= group_num+1 & i <= 2*group_num){
          tsne_plot$x$data[[i]]$text <- paste("\nGroups :",
                                              tsne_plot$x$data[[i]]$name)
        }
      }
      for (i in seq_len((2*group_num))) {
        rev_n <- rev(seq_len((2*group_num)))
        if(i == 1){rev_list <- list(tsne_plot$x$data[[rev_n[i]]])}
        else{rev_list<-c(rev_list, list(tsne_plot$x$data[[rev_n[i]]]))}
      }
      tsne_plot$x$data <- rev_list
    }

  }

  return(list(tsne_data, tsne_plot))
}




