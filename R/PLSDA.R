#' @title PLSDA
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
#' @param ncomp A numeric value. The number of components to include in 
#' the model. (default: 2)
#' @param scaling A logical value. If scaling = TRUE, each block is 
#' standardized to zero means and unit variances. (default: TRUE).
#' @param cluster_method A character string indicating which method to be used
#'  for clustering. Allowed method include \bold{"kmeans"}, \bold{"kmedoids"}, 
#'  \bold{"hclustering"}, \bold{"dbscan"}, \bold{"group_info"}.
#' @param group_num A positive integer specifying the number of clusters. 
#' Group number must be between 1 and 10.
#' @param var1 A character string or numeric. According to the 
#' selected cluster method.
#' \itemize{
#'   \item \emph{kmeans}: A specifying the metric to be used for calculating 
#'   dissimilarities between observations include \bold{"euclidean"} 
#'   and \bold{"manhattan"}.
#'   \item \emph{hclustering}: which correlation coefficient is to be computed 
#'   or which the distance measure to be used include \bold{"pearson"}, 
#'   \bold{"kendall"}, \bold{"spearman"}, \bold{"euclidean"}, 
#'   \bold{"manhattan"}, \bold{"maximum"}, \bold{"canberra"}, 
#'   \bold{"binary"}, \bold{"minkowski"}.
#'   \item \emph{dbscan}: size of the epsilon neighborhood.
#' }
#' @param var2 A character string or numeric. According to the 
#' selected cluster method.
#' \itemize{
#'   \item \emph{hclustering}: The agglomeration method to be used. 
#'   This should be (an unambiguous abbreviation of) one of \bold{"ward.D"}, 
#'   \bold{"ward.D2"}, \bold{"single"}, \bold{"complete"}, 
#'   \bold{"average"}(= UPGMA), \bold{"mcquitty"}(= WPGMA), 
#'   \bold{"median"}(= WPGMC), or \bold{"centroid"}(= UPGMC).
#'   \item \emph{dbscan}: number of minimum points in the 
#'   eps region (for core points).
#' }
#' @param insert_ref_group A character string. The name of 'ctrl' after 
#' name conversion.
#' @param ref_group A character string. The name of 'exp' after 
#' name conversion.
#' @return Return a list with 2 data frames and 2 plots.
#' \enumerate{
#' \item X.variate: a data frame of sample variate associated to each 
#' component if PLS-DA is applied. This is for \bold{sample.plot}.
#' \item X.loading: a data frame of sample loading, which are coefficients 
#' assigned to each variable to define each component. These coefficients 
#' indicate the importance of each variable in PLS-DA. This plot is for 
#' \bold{variable.plot}.
#' \item sample.plot: PLS-DA sample plot. The axis labels indicate the 
#' amount of variation explained per component.
#' \item variable.plot: PLS-DA variable plot, display the variables that 
#' contribute to the definition of each component.  These variables should be 
#' located towards the circle of radius 1, far from the center.
#' }
#' @export
#' @examples
#' data("DE_exp_data")
#' data("DE_group_info")
#' exp_data <- DE_exp_data
#' group_info <- DE_group_info
#' exp_transform_table <- data_process(exp_data, exclude_var_missing=TRUE,
#'                                     missing_pct_limit=50, replace_zero=TRUE,
#'                                     zero2what='min', xmin=0.5,
#'                                     replace_NA=TRUE, NA2what='min',
#'                                     ymin=0.5,  pct_transform=TRUE,
#'                                     data_transform=FALSE, trans_type='log',
#'                                     centering=FALSE, scaling=FALSE)
#' DE_species_table_sig <- DE_species_2(exp_transform_table,
#'                                      data_transform=TRUE,
#'                                      group_info=group_info, paired=FALSE,
#'                                      test='t.test', adjust_p_method='BH',
#'                                      sig_stat='p.adj', sig_pvalue=0.05,
#'                                      sig_FC=2)$DE_species_table_sig
#' PLSDA(exp_transform_table, group_info=group_info,
#'       sig_feature=DE_species_table_sig$feature,
#'       ncomp=2, scaling=TRUE, cluster_method='group_info',
#'       group_num=NULL, var1=NULL, var2=NULL,
#'       insert_ref_group=NULL, ref_group=NULL)
PLSDA <- function(exp_transform_table, group_info=NULL, sig_feature=NULL, 
                  ncomp=2, scaling=TRUE, cluster_method, group_num=NULL, 
                  var1=NULL, var2=NULL, insert_ref_group=NULL, 
                  ref_group=NULL, plotly=TRUE){
  
  exp_transform_table <- as.data.frame(exp_transform_table)
  if(!is.null(group_info)){
    group_info <- as.data.frame(group_info)
  }
  if(ncol(exp_transform_table) < 7){
    stop("At least 6 samples.")
  }
  if(nrow(exp_transform_table) < 7){
    stop("At least 6 features.")
  }
  if(sum(vapply(exp_transform_table[, -1], class, character(1)) %in% 
         c("numeric", "integer")) != ncol(exp_transform_table[, -1])){
    stop("First column type must be 'character',others must be 'numeric'")
  }
  if(tibble::is.tibble(exp_transform_table)){
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
      if(sum(group_info[, 1] %in% colnames(exp_transform_table)) != 
         nrow(group_info) | 
         sum(group_info[, 1] %in% colnames(exp_transform_table)) != 
         ncol(exp_transform_table[, -1])){
        stop("group_info 'sample_name' must same as the name of samples of 
             exp_transform_table")
      }
      if(length(unique(group_info[, 3])) == 2){
        if(sum(table(group_info[, 3]) >= 1) != 2){
          stop("group_info column 'group' only can have 2 groups, and >= 1 
               sample for each group.")
        }
      }else{
        stop("group_info column 'group' only can have 2 groups, and >= 1 
             sample for each group.")
      }
      if(!is.null(insert_ref_group)){
        if(!insert_ref_group %in% group_info[,3]){
          stop("The reference group entered by users must be included in the group_info.")
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
      if(sum(group_info[, 1] %in% colnames(exp_transform_table)) != 
         nrow(group_info) | 
         sum(group_info[, 1] %in% colnames(exp_transform_table)) != 
         ncol(exp_transform_table[, -1])){
        stop("group_info 'sample_name' must same as the name of samples of 
             exp_transform_table")
      }
      if(length(unique(group_info[, 3])) == 2){
        if(sum(table(group_info[, 3]) >= 1) != 2){
          stop("group_info column 'group' only can have 2 groups, and >= 1 
               sample for each group.")
        }
      }else{
        stop("group_info column 'group' only can have 2 groups, and >= 1 
             sample for each group.")
      }
      if(!is.null(insert_ref_group)){
        if(!insert_ref_group %in% group_info[, 3]){
          stop("The insert_ref_group entered by users must be included in 
               the group_info.")
        }
      }
    }
  }

  if(cluster_method=='kmeans'){
    if(!is.null(group_num)){
      if( group_num < 1 ||  group_num > 10){
        stop('Group number must be between 1 and 10')
      }
    }else{
      stop('Group number can not be NULL')
    }
  }else if(cluster_method == 'kmedoids'){
    if(!var1 %in% c("euclidean", "manhattan")){
      stop('var1 must be one of the strings "euclidean" or "manhattan".')
    }
  }else if(cluster_method == 'hclustering'){
    if(!var1 %in% c('pearson', 'spearman', 'kendall', "euclidean", "maximum", 
                    "manhattan", "canberra", "binary", "minkowski")){
      stop('var1 must be one of the strings "pearson", "kendall", "spearman", 
           "euclidean", "maximum", "manhattan", "canberra", 
           "binary", or "minkowski"')
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
  if(!is.null(insert_ref_group) & !is.null(ref_group) &!is.null(group_info)){
    exp_raw_name <- ref_group[-which(insert_ref_group %in% ref_group)]
    group_info$group[which(group_info$group == 'ctrl')] <-  insert_ref_group
    group_info$group[which(group_info$group == 'exp')] <-  exp_raw_name
  }

  colnames(exp_transform_table)[1] <- 'feature'

  if(sum(!colnames(exp_transform_table) %in% group_info$sample_name) != 1){
    stop("Sample name not matched in the exp_transform_table and group_info!")
  }

  if(!is.null(sig_feature)){
    exp_transform_table <- exp_transform_table %>%
      dplyr::filter(eval(parse(
        text=colnames(exp_transform_table)[1])) %in% sig_feature)
  }

  num <- apply(exp_transform_table[-1], 1, FUN=function(x){length(unique(x))})
  exp_transform_table <- exp_transform_table[(num!=1), ]
  exp_transform_table <- exp_transform_table[!is.infinite(
    rowSums(exp_transform_table[-1], na.rm=TRUE)), ]

  group_info <- group_info %>% 
    dplyr::filter(sample_name %in% colnames(exp_transform_table)[-1])

  if(!is.null(insert_ref_group) & !is.null(ref_group) &!is.null(group_info)){
    Y <- ifelse(group_info$group == insert_ref_group, 1, 2)
  }else{
    Y <- ifelse(group_info$group == 'ctrl', 1, 2)
  }

  names(Y) <- group_info$sample_name
  rownames(exp_transform_table) <- exp_transform_table$feature

  X <- exp_transform_table[c('feature', group_info$sample_name)] %>%
    dplyr::select(-1) %>%
    t()

  plsda.res <- mixOmics::plsda(X, Y, ncomp, scaling)
  X.variate <- plsda.res$variates$X
  X.loading <- plsda.res$loadings$X

  #grouping
  if(cluster_method == 'kmeans'){
    cluster_group <- stats::kmeans(X.variate, centers=group_num)$cluster
  }else if(cluster_method == 'kmedoids'){
    cluster_group <- cluster::pam(X.variate, k=group_num, 
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
    cluster_group <- stats::hclust(dist.fun(X.variate), method=var2)
    cluster_group <- stats::cutree(cluster_group, k=group_num)
  }else if(cluster_method == 'dbscan'){
    cluster_group <- dbscan::dbscan(X.variate, eps=var1, minPts=var2)$cluster
    cluster_group <- ifelse(cluster_group>9,9,cluster_group) %>% as.character()
    cluster_group <- ifelse(cluster_group == '0', 'noise', cluster_group)
  }else if(cluster_method == 'group_info'){
    group_order <- purrr::map_dbl(rownames(X.variate), 
                                  function(x){which(
                                    x == group_info$sample_name)})

    cluster_group <- group_info[group_order,] %>% .$group
    group_num <- length(unique(cluster_group))

  }

  cluster_group <- as.character(cluster_group)

  color <- c("#00AFBB", "#E7B800", "#FC4E07", "#42B540FF", "#BB3099", 
             "#EE0099", "#0000AC", "#868686FF", '#00468BFF', 'black')
  if(cluster_method %in% c('dbscan')){
    sample.plot <-  factoextra::fviz_cluster(
      list(data=as.data.frame(X.variate), cluster=cluster_group), 
      palette=color[seq_len(length(unique(cluster_group)))], 
      geom="point", stand=FALSE, ellipse=FALSE, ggtheme=ggplot2::theme_bw(), 
      show.clust.cent=FALSE, xlab='PLSDA-1', ylab='PLSDA-2', main="PLSDA",
      legend.title='Groups')
    if(plotly == TRUE){
      sample_plot <- plotly::ggplotly(sample.plot)
      for (i in seq_len(length(sample_plot$x$data))){
        if (!is.null(sample_plot$x$data[[i]]$name)){
          sample_plot$x$data[[i]]$name <- gsub(
            "\\(", "", 
            stringr::str_split(sample_plot$x$data[[i]]$name, ",")[[1]][1])
        }
        if(i<=length(unique(cluster_group))){
          sample_plot$x$data[[i]]$text <- paste(
            "x :", round(sample_plot$x$data[[i]]$x, 3),
            "\ny :", round(sample_plot$x$data[[i]]$y, 3),
            "\nGroups :", sample_plot$x$data[[i]]$name,
            "\nSample name :", sample.plot$data[which(
              sample.plot$data$cluste == sample_plot$x$data[[i]]$name), ]$name)
        }
      }
      sample.plot <- sample_plot
    }
    
  }else{
    sample.plot <- factoextra::fviz_cluster(
      list(data=as.data.frame(X.variate), 
           cluster=cluster_group), palette=color[seq_len(group_num)], 
      geom="point", stand=FALSE, ellipse=TRUE, ellipse.type='norm', 
      ggtheme=ggplot2::theme_bw(), show.clust.cent=FALSE, xlab='PLSDA-1', 
      ylab='PLSDA-2', main="PLSDA", legend.title='Groups')
    if(plotly == TRUE){
      sample_plot <- plotly::ggplotly(sample.plot)
      for (i in seq_len(length(sample_plot$x$data))){
        if (!is.null(sample_plot$x$data[[i]]$name)){
          sample_plot$x$data[[i]]$name <- gsub(
            "\\(", "", 
            stringr::str_split(sample_plot$x$data[[i]]$name, ",")[[1]][1])
        }
        if (i <= group_num){
          sample_plot$x$data[[i]]$text <- paste(
            "x :", round(sample_plot$x$data[[i]]$x, 3), 
            "\ny :",round(sample_plot$x$data[[i]]$y, 3), 
            "\nGroups :", sample_plot$x$data[[i]]$name, 
            "\nSample name :", sample.plot$data[which(
              sample.plot$data$cluster == sample_plot$x$data[[i]]$name), ]$name)
        }else if(i>=group_num+1 & i <= 2*group_num){
          sample_plot$x$data[[i]]$text <- paste("Groups :", 
                                                sample_plot$x$data[[i]]$name)
        }
      }
      for (i in seq_len((2*group_num))) {
        rev_n <- rev(seq_len((2*group_num)))
        if(i == 1){
          rev_list <- list(sample_plot$x$data[[rev_n[i]]])
        }
        else{rev_list<-c(rev_list, list(sample_plot$x$data[[rev_n[i]]]))}
      }
      sample_plot$x$data <- rev_list
      sample.plot <- sample_plot
    }
  }


  rown <- rownames(X.variate)

  X.variate <- X.variate %>% as.data.frame() %>%
    dplyr::mutate(sample_name=rown) %>%
    dplyr::mutate(group=cluster_group) %>%
    dplyr::select(sample_name, group, dplyr::everything())
  colnames(X.variate)[3:4] <- c('PLSDA1', 'PLSDA2')

  X.loading <- X.loading %>% as.data.frame()
  colnames(X.loading) <- c('PLSDA1', 'PLSDA2')

  variable.tab <- mixOmics::plotVar(plsda.res, comp=seq_len(2), 
                                    var.names=TRUE, plot=FALSE)
  if(plotly == TRUE){
    creat_circle_data <- function(r, x=0, y=0){
      circle.x <- seq(x-r, x+r, length.out =150)
      circle.y <- c(sqrt(r^2-circle.x^2), -sqrt(r^2-circle.x^2))
      circle.x <- c(circle.x,seq(x+r, x-r, length.out =150))
      circle <- data.frame(circle.x, circle.y)
    }
    circle_data_1 <- creat_circle_data(1)
    circle_data_2 <- creat_circle_data(0.5)
    
    variable.plot <- plotly::plot_ly(x=~circle.x, 
                                     y=~circle.y, 
                                     hoverinfo="none") %>%
      plotly::add_trace(data=circle_data_1, 
                        type="scatter", 
                        mode="lines", hoverinfo="none", showlegend=FALSE, 
                        ine=list(color='black', width=2)) %>%
      plotly::add_trace(data=circle_data_2, 
                        type="scatter", mode="lines", 
                        hoverinfo="none", showlegend=FALSE, 
                        line=list(color='black', width=2)) %>%
      plotly::add_trace(data=variable.tab, x=~x, y=~y, type="scatter", 
                        mode="markers", hoverinfo="text", showlegend=FALSE, 
                        marker=list(color="red", size=5), 
                        text=~paste("x :", round(x, 3), 
                                    "\ny :", round(y, 3), 
                                    "\nfeature :", variable.tab$names)) %>%
      plotly::layout(xaxis=list(title="PLSDA-1"), 
                     yaxis=list(title="PLSDA-2"), 
                     title='Variables - PLSDA')
  }else{
    circles <- data.frame(x0 = rep(0, 2),
                          y0 = rep(0, 2),
                          r = c(0.5, 1), 
                          stringsAsFactors = F)
    
    variable.plot <- ggplot() +
      geom_circle(data=circles, 
                  aes(x0=x0, y0=y0, r=r), show.legend = F) + 
      geom_text(data=variable.tab, 
                aes(x=x, y=y, label=names), colour='#8B0000', 
                check_overlap=F, size=4, show.legend=F) + 
      theme_bw() + 
      labs(x='PLSDA-1', y='PLSDA-2', title='Variables - PLSDA') + 
      theme(plot.title=element_text(hjust=0.5))
  }
  return(list(X.variate, X.loading, sample.plot, variable.plot))

}
