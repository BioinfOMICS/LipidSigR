#' @title PCA
#' @description PCA (Principal Component Analysis):  A method of dimensionality
#' reduction to transform data from a high-dimensional space into a
#' low-dimensional space so that the vital properties of the original data
#' can be retained.
#' @description This function calculates Principal component analysis with
#' the classical prcomp function and visualizes results by PCA plot and scree
#' plot. Furthermore, this function also visualizes the features' contribution
#' of the user-defined principal components and the correlation between a
#' feature (lipid species) and a principal component (PC).
#' @param exp_transform_table A data frame of predictors, including
#' features (molecule, lipid class, etc.) and their expression of each sample.
#' NAs are not allowed. The name of the first column must be
#' "feature" (lipid species).
#' @param group_info A data frame with sample name, label name, the group name
#' of samples (ctrl & exp), and pair number representing 'the pair' for the
#' t-test/Wilcoxon test. NAs are allowed.
#' @param sig_feature A character string for identifying significant variable
#' in \bold{exp_transform_table}.
#' @param scaling A logical value. If scaling = TRUE, each block is
#' standardized to zero means and unit variances. (default: TRUE).
#' @param centering A logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal the
#' number of columns of x can be supplied. The value is passed to scale.
#' @param cluster_method A character string indicating which method to be used
#' for clustering. Allowed method include \bold{'kmeans'}, \bold{'kmedoids'},
#' \bold{'hclustering'}, \bold{'dbscan'},
#' \bold{'group_info'}. (default: 'kmeans')
#' @param group_num A positive integer specifying the number of clusters.
#' The number must be between 1 and 10.
#' @param var1 A character string or numeric. According to the selected
#' cluster method.
#' \itemize{
#'   \item \emph{kmeans}: specifying the metric to be used for calculating
#'   dissimilarities between observations include \bold{"euclidean"} and
#'   \bold{"manhattan"}.
#'   \item \emph{hclustering}: which correlation coefficient is to be computed
#'   or which the distance measure to be used include \bold{"pearson"},
#'   \bold{"kendall"}, \bold{"spearman"}, \bold{"euclidean"},
#'   \bold{"manhattan"}, \bold{"maximum"}, \bold{"canberra"},
#'   \bold{"binary"}, \bold{"minkowski"}.
#'   \item \emph{dbscan}: size of the epsilon neighborhood.
#' }
#' @param var2 A character string or numeric. According to the selected
#' cluster method.
#' \itemize{
#'   \item \emph{hclustering}: The agglomeration method to be used. This
#'   should be (an unambiguous abbreviation of) one of \bold{"ward.D"},
#'   \bold{"ward.D2"}, \bold{"single"}, \bold{"complete"},
#'   \bold{"average"} (=UPGMA), \bold{"mcquitty"} (= WPGMA),
#'   \bold{"median"} (= WPGMC), or \bold{"centroid"} (= UPGMC).
#'   \item \emph{dbscan}: number of minimum points in the eps
#'   region (for core points).
#' }
#' @param insert_ref_group A character string. The name of 'ctrl' after
#' name conversion.
#' @param ref_group A character string. The name of 'exp' after name conversion.
#' @param n_PC A numeric vector specifying the dimension(s) of interest.
#' @param top_n_feature A numeric value specifying the number of top
#' elements to be shown.
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @return Return a list with 1 model, 2 data frame, and 4 plots.
#' \enumerate{
#' \item pca: prcomp (lists of the standard deviations, the matrix of variable
#' loading, the value of the rotated data, the centering and scaling )
#' \item pca_rotated_data: a data frame of PCA rotated data
#' \item pca_contrib_table: data frame, PCA contribution table
#' \item pca_biplot: PCA plot
#' \item pca_screeplot: Scree plot of top n principle components
#' \item feature_contrib: plot, the contribution of top N features of the
#' user-defined principal components.
#' \item pca_variable: correlation circle plot(factor map) of PCA variables.
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
#' PCA(exp_transform_table, group_info=NULL, sig_feature=NULL, scaling=TRUE,
#'     centering=TRUE, cluster_method='kmeans', group_num=2, var1=NULL,
#'     var2=NULL, insert_ref_group=NULL,
#'     ref_group=NULL, n_PC=c(1,2), top_n_feature=10, plotly=TRUE)
PCA <- function(exp_transform_table, group_info=NULL, sig_feature=NULL,
                scaling=TRUE, centering=TRUE, cluster_method='kmeans',
                group_num=NULL, var1=NULL, var2=NULL,
                insert_ref_group=NULL, ref_group=NULL,
                n_PC, top_n_feature, plotly=TRUE){

  exp_transform_table <- as.data.frame(exp_transform_table)
  if(!is.null(group_info)){
    group_info <- as.data.frame(group_info)
  }
  if(ncol(exp_transform_table) == 2){
    if(sum(class(exp_transform_table[, -1]) %in% c("numeric", "integer")) != 1){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(sum(vapply(exp_transform_table[, -1], class, character(1)) %in%
           c("numeric", "integer")) != ncol(exp_transform_table[, -1])){
      stop("First column type must be 'character',others must be 'numeric'")
    }
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
  if(ncol(exp_transform_table) < 7){
    stop("At least 6 samples.")
  }
  if(nrow(exp_transform_table) < 3){
    stop("At least 2 features.")
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
      if(sum(vapply(group_info[, seq_len(3)],
                    class,character(1)) != "character") == 0){
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
      if(sum(!is.na(group_info[, 4])) !=0 |
         sum(table(group_info[, 4]) != 2) != 0 &
         sum(is.na(group_info[, 4])) != 0){
        stop("group_info each pair must have a specific number, staring from
             1 to N. Cannot have NA, blank, or skip numbers.")
      }
      if(sum(group_info[, 1] %in% colnames(exp_transform_table)) !=
         nrow(group_info) | sum(group_info[, 1] %in%
                                colnames(exp_transform_table)) !=
         ncol(exp_transform_table[, -1])){
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
      if(sum(vapply(group_info, class, character(1)) != "character") != 0){
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
        if(sum(table(group_info[, 3])>=1) != 2){
          stop("group_info column 'group' only can have 2 groups, and >= 1
               sample for each group.")
        }
      }else{
        stop("group_info column 'group' only can have 2 groups, and >= 1 sample
             for each group.")
      }
      if(!is.null(insert_ref_group)){
        if(!insert_ref_group %in% group_info[, 3]){
          stop("The insert_ref_group entered by users must be included in the
               group_info.")
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
      dplyr::filter(eval(
        parse(text=colnames(exp_transform_table)[1])) %in% sig_feature)
  }

  num <- apply(exp_transform_table[-1], 1, FUN=function(x){length(unique(x))})
  exp_transform_table <- exp_transform_table[(num!=1),]
  exp_transform_table <- exp_transform_table[!is.infinite(rowSums(
    exp_transform_table[-1], na.rm=TRUE)),]

  pca_table <- exp_transform_table[-1] %>% t() %>% as.data.frame()

  colnames(pca_table) <- exp_transform_table[[1]]


  pca_table <- stats::na.omit(pca_table)
  pca <- stats::prcomp(pca_table, scale=scaling, center=centering)

  #output pca data
  pca_rotated_data <- as.data.frame(pca$x) %>%
    dplyr::mutate(sample_name=colnames(exp_transform_table)[-1])


  #grouping
  if(cluster_method == 'kmeans'){
    cluster_group <- stats::kmeans(pca$x, centers=group_num)$cluster
  }else if(cluster_method == 'kmedoids'){
    cluster_group <- cluster::pam(pca$x, k=group_num,
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
    cluster_group <- stats::hclust(dist.fun(pca$x), method=var2)
    cluster_group <- stats::cutree(cluster_group, k=group_num)
  }else if(cluster_method == 'dbscan'){
    cluster_group <- dbscan::dbscan(pca$x, eps=var1, minPts=var2)$cluster
    cluster_group <- ifelse(cluster_group > 9, 9,
                            cluster_group) %>% as.character()
    cluster_group <- ifelse(cluster_group == '0', 'noise', cluster_group)
  }else if(cluster_method == 'group_info'){

    group_order <- purrr::map_dbl(rownames(pca_table),
                                  function(x){which(
                                    x == group_info$sample_name)})

    cluster_group <- group_info[group_order, ] %>% .$group
    group_num <- length(unique(cluster_group))
  }


  pca_rotated_data <- pca_rotated_data %>%
    dplyr::mutate(group=cluster_group) %>%
    dplyr::select(sample_name, group, dplyr::everything())
  if(!is.null(insert_ref_group) & !is.null(ref_group)){
    exp_raw_name <- ref_group[-which(insert_ref_group %in% ref_group)]
    pca_rotated_data$group[which(pca_rotated_data$group == 'ctrl')] <-
      insert_ref_group
    pca_rotated_data$group[which(pca_rotated_data$group == 'exp')] <-
      exp_raw_name
  }
  if(TRUE){


    #pca contribution table
    pca_contrib_table <- factoextra::get_pca_var(pca)$contrib %>%
      as.data.frame() %>%
      dplyr::mutate(feature=rownames(.)) %>%
      dplyr::select(feature, dplyr::everything())

    PC <- ncol(pca_contrib_table)-1
    colnames(pca_contrib_table)[seq(2, ncol(pca_contrib_table))] <-
      stringr::str_c('PC', as.character(seq_len(PC)))

    #pca plot*5
    color <- c("#00AFBB", "#E7B800", "#FC4E07", "#42B540FF", "#BB3099",
               "#EE0099", "#0000AC", "#868686FF", '#00468BFF', 'black')

    if(cluster_method %in% c('dbscan')){
      pca_biplot <- factoextra::fviz_pca_ind(
        pca,label="none",
        habillage=cluster_group,
        palette=color[seq_len(length(unique(cluster_group)))],
        addEllipses=FALSE, invisible="quali",
        title='PCA')
      if(plotly == TRUE){
        biplot <- plotly::ggplotly(pca_biplot)
        for (i in seq_len(length(biplot$x$data))){
          if (!is.null(biplot$x$data[[i]]$name)){
            biplot$x$data[[i]]$name <- gsub(
              "\\(", "",
              stringr::str_split(biplot$x$data[[i]]$name, ",")[[1]][1])
          }
          if(i <= length(unique(cluster_group))){
            biplot$x$data[[i]]$text <- paste(
              "x :",  round(biplot$x$data[[i]]$x, 3),
              "\ny :", round(biplot$x$data[[i]]$y, 3),
              "\nGroups :", biplot$x$data[[i]]$name,
              "\nSample name :", pca_biplot$data[which(
                pca_biplot$data$Groups == biplot$x$data[[i]]$name),]$name)
          }
        }
        pca_biplot <- biplot
      }

    }else{
      pca_biplot <- factoextra::fviz_pca_ind(pca, label="none",
                                             habillage=cluster_group,
                                             palette=color[seq_len(group_num)],
                                             addEllipses=TRUE,
                                             invisible="quali",
                                             title='PCA')
      if(plotly == TRUE){
        biplot <- plotly::ggplotly(pca_biplot)
        for (i in seq_len(length(biplot$x$data))){
          if (!is.null(biplot$x$data[[i]]$name)){
            biplot$x$data[[i]]$name <- gsub(
              "\\(", "",
              stringr::str_split(biplot$x$data[[i]]$name, ",")[[1]][1])
          }
          if (i <= group_num){
            biplot$x$data[[i]]$text <- paste("Dim1 :",
                                             round(biplot$x$data[[i]]$x, 3),
                                             "\nDim2 :",
                                             round(biplot$x$data[[i]]$y, 3),
                                             "\nGroups :",
                                             biplot$x$data[[i]]$name,
                                             "\nSample name :",
                                             pca_biplot$data[which(
                                               pca_biplot$data$Groups ==
                                                 biplot$x$data[[i]]$name)
                                               ,]$name)
          }else if(i >= group_num+1 & i <= 2*group_num){
            biplot$x$data[[i]]$text <- paste("\nGroups :",
                                             biplot$x$data[[i]]$name)
          }
        }
        for (i in seq_len((2*group_num))) {
          rev_n <- rev(seq_len((2*group_num)))
          if( i== 1){rev_list <- list(biplot$x$data[[rev_n[i]]])}
          else{rev_list <- c(rev_list, list(biplot$x$data[[rev_n[i]]]))}
        }
        biplot$x$data <- rev_list
        pca_biplot <- biplot
      }
    }

    pca_screeplot <- factoextra::fviz_screeplot(
      pca,
      ncp=10,
      main="Scree plot for Top10 PCs",
      xlab="Principle components (PCs)",
      ylab='Explained variance (%)')
    if(plotly == TRUE){
      pca_screeplot <- plotly::ggplotly(pca_screeplot) %>%
        plotly::style(traces=seq_len(2), hoverinfo="text",
                      text=paste("Principle components :",
                                 pca_screeplot$data$dim,
                                 "\nExplained variance :",
                                 round(pca_screeplot$data$eig, 3),"%"))
    }

  }
  feature_contrib <- factoextra::fviz_contrib(
    pca, choice="var", addlabels=TRUE, axes=n_PC, top=top_n_feature,
    title=stringr::str_c('Contribution of Top', as.character(top_n_feature),
                         ' features to PC-', stringr::str_c(as.character(n_PC),
                                                            collapse=',')))
  if(plotly == TRUE){
    feature_contrib <- plotly::ggplotly(feature_contrib)
    feature_contrib$x$data[[1]]$text <- paste(
      "Feature :", feature_contrib$x$layout$xaxis$categoryarray,
      "\nContributions :", round(feature_contrib$x$data[[1]]$y, 3),"%")
    feature_contrib$x$data[[2]]$text <- paste(
      "yintercept :", round(feature_contrib$x$data[[2]]$y[1], 3))
  }

  eig.val <- factoextra::get_eigenvalue(pca)

  pca_variable <- factoextra::fviz_pca_var(
    pca, col.var="contrib", select.var=list(contrib=top_n_feature),
    gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"),
    repel=TRUE, title='Variables - PCA',
    xlab=stringr::str_c("PC1 (", as.character(round(eig.val[1, 2], 1)), '%)'),
    ylab=stringr::str_c("PC2 (", as.character(round(eig.val[2, 2], 1)), '%)'))
  if(plotly == TRUE){
    pca_variable$data <- dplyr::arrange(pca_variable$data,contrib)
    pca_variable_ggplotly <- plotly::ggplotly(pca_variable)
    for (i in seq_len(top_n_feature)) {
      pca_variable_ggplotly <- pca_variable_ggplotly %>%
        plotly::add_annotations(
          data=pca_variable$data[i,],
          x=~x, y=~y, text="",
          arrowcolor=pca_variable_ggplotly$x$data[[i+1]]$line$color,
          showarrow=TRUE, axref='x', ayref='y', ax=0, ay=0)
      pca_variable_ggplotly$x$data[[i+1]]$text <-
        paste("X :", round(pca_variable$data$x[i], 3),
              "\nY :", round(pca_variable$data$y[i], 3),
              "\nName :", pca_variable$data$name[i])
    }
    pca_variable_ggplotly$x$data[[top_n_feature+2]]$hoverinfo <- "none"
    pca_variable <-  pca_variable_ggplotly
  }
  if(top_n_feature > ncol(pca$rotation)){
    top_n_feature <- ncol(pca$rotation)
    warning_message <- paste0("top_n_feature will only show ",
                              ncol(pca$rotation))
    warning(warning_message)
  }


  return(list(pca, pca_rotated_data, pca_contrib_table,
              pca_biplot, pca_screeplot,feature_contrib,pca_variable))

}



