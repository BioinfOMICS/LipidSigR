#' @title corr_heatmap
#' @description The correlation heatmap illustrates the correlation between
#' lipid species or samples and depicts the patterns in each group.
#' The correlation can be calculated by Pearson or Spearman.
#' \enumerate{
#' \item Compute and cluster the correlation coefficient depending on
#' user-defined method and distance.
#' \item Print two correlation heatmaps to illustrate the correlation between
#' lipid species or samples.
#' }
#' @param exp_data A data frame of predictors, including features
#' (molecule, lipid class, etc.) and their expression of each sample. NAs are
#' not allowed. The name of the first column must be "feature" (lipid species).
#' @param corr_method A character string indicating which correlation
#' coefficient is to be computed. One of "pearson" or
#' "spearman". (default: "pearson")
#' @param distfun A character string of the distance measure indicating which
#' correlation coefficient (or covariance) is to be computed. Allowed methods
#' include "pearson", "spearman", "kendall", euclidean", "maximum",
#' "manhattan", "canberra", "binary", "minkowski". (default: "maximum")
#' @param hclustfun A character string of the agglomeration method to be used.
#' This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). (default: "average")
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @param type It should be 'sample' or 'lipid.' 'sample' outputs the
#' correlation results of samples, and 'lipid' outputs output the correlation
#' results of lipid species.
#' @return Return 1 plot and a list of 3 matrices.
#' \enumerate{
#' \item corr_coef: a matrix of correlation coefficients between samples
#' or lipids.
#' \item corr_p: a matrix of correlation p-value between samples or lipids.
#' \item heatmap: a heatmap plot by samples or lipid species.
#' \item reorder_data: a reorder sample or lipids correlation coefficients
#' matrix.
#' }
#' @export
#' @examples
#' data("profiling_exp_data")
#' exp_data <- profiling_exp_data[1:50, ]
#' exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
#'                               missing_pct_limit=50,
#'                               replace_zero=TRUE, zero2what='min',
#'                               xmin=0.5, replace_NA=TRUE,
#'                               NA2what='min', ymin=0.5, pct_transform=TRUE,
#'                               data_transform=TRUE, trans_type='log',
#'                               centering=FALSE, scaling=FALSE)
#' corr_heatmap(exp_transform, corr_method="pearson", distfun="maximum",
#'              hclustfun="average", plotly=TRUE, type='sample')
corr_heatmap <- function(exp_data,corr_method="pearson",
                         distfun="maximum", hclustfun='average',
                         plotly=TRUE,type='sample'){
  exp_data <- as.data.frame(exp_data)
  if(!is(exp_data[, 1], 'character')){
    stop("The first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data) == 2){
    if(!is(exp_data[, 1], 'character') |
       sum(class(exp_data[, -1]) %in% c("numeric", "integer")) != 1){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(!is(exp_data[, 1], 'character') |
       sum(vapply(exp_data[, -1], class, character(1)) %in%
           c("numeric", "integer")) != ncol(exp_data[, -1])){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }
  if(tibble::is_tibble(exp_data)){
    if(nrow(exp_data) != nrow(unique(exp_data[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(exp_data) != length(unique(exp_data[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if(ncol(exp_data) < 3){
    stop("At least 2 samples.")
  }
  if(sum(is.na(exp_data[, -1])) > 0){
    stop("Variables can not be NA")
  }
  
  cor.test.p <- function(x,method){
    FUN <- function(x, y,method) stats::cor.test(x=x,
                                                 y=y,
                                                 method=method,
                                                 na.action=
                                                   "na.exclude")[["p.value"]]
    z <- outer(
      colnames(x),
      colnames(x),
      Vectorize(function(i, j) FUN(x[, i], x[, j], method=method))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
  }
  
  if(type == 'sample'){
    exp_data <- subset(exp_data, select=-feature)
    corr_coef <- stats::cor(exp_data, method=corr_method,
                            use="pairwise.complete.obs")
    corr_p <- cor.test.p(x=exp_data, method=corr_method)
    data <- corr_coef
  }else if(type == 'lipid'){
    exp_data <- exp_data %>%
      tidyr::gather(-feature, key='sample_name', value='value') %>%
      tidyr::spread(key='feature', value='value')
    exp_data <- subset(exp_data, select=-sample_name)
    corr_coef <- stats::cor(exp_data, method=corr_method,
                            use="pairwise.complete.obs")
    corr_p <- cor.test.p(x=exp_data, method=corr_method)
    data <- corr_coef
  }
  
  
  if(sum(is.na(data)) == 0 & nrow(data) >= 2 & ncol(data) >= 2){
    if(plotly == TRUE){
      if(distfun %in% c('pearson','spearman','kendall')){
        dist_fun <- function(x){
          x <- t(x)
          cor.mat <- stats::cor(x, method=distfun, use='complete.obs')
          cor.mat <- (1-cor.mat)
          cor.dist <- stats::as.dist(cor.mat)
          return(cor.dist)
        }
      }else{
        dist_fun <- function(x) stats::dist(x, method=distfun)
      }
      hclust_fun <- function(x) stats::hclust(x, method=hclustfun) 
      if(min(data) >= 0 || max(data) <= 0){
        heatmap <- heatmaply::heatmaply(
          data, distfun=dist_fun, hclustfun=hclust_fun,
          scale_fill_gradient_fun=
            ggplot2::scale_fill_gradient2(low="#0571b0", mid="white", 
                                          high="#ca0020", midpoint=0, 
                                          na.value="grey50"),
          column_text_angle=90, grid_color="white", 
          margins=c(l=0.2, r=0, t=0, b=80),
          hide_colorbar=T)
      }else{
        heatmap <- heatmaply::heatmaply(
          data, distfun=dist_fun, hclustfun=hclust_fun,
          scale_fill_gradient_fun=
            ggplot2::scale_fill_gradient2(low="#0571b0", mid="white", 
                                          high="#ca0020", midpoint=0, 
                                          na.value="grey50"),
          column_text_angle=90, grid_color="white", 
          margins=c(l=0.2, r=0, t=0, b=80))
      }
      
      ax <- list(title="", zeroline=FALSE, showline=FALSE, 
                 showgrid=FALSE)
      bx <- list(title="", zeroline=FALSE, showline=FALSE, 
                 showticklabels=FALSE, showgrid=FALSE, ticks='')
      
      if(nrow(data) > 50){
        heatmap <- heatmap  %>%
          layout(xaxis=bx) 
      }else{
        heatmap <- heatmap  %>%
          layout(xaxis=ax) 
      }
      if(ncol(data) > 50){
        heatmap <- heatmap  %>%
          layout(yaxis2=bx) 
      }else{
        heatmap <- heatmap  %>%
          layout(yaxis2=ax) 
      }
      ## reorder Sample correlation matrix ##
      reorder_data <- data[rev(heatmap$x$layout$yaxis2$ticktext), 
                           heatmap$x$layout$yaxis2$ticktext]
    }else{
      if(distfun %in% c("pearson", "kendall", "spearman")){
        col_dend <- stats::hclust(
          stats::as.dist(1-stats::cor(data, method=distfun)),
          method=hclustfun)
      }else{
        col_dend <- stats::hclust(
          stats::dist(t(data), method=distfun), method=hclustfun)
      }
      sample_col_dend <- function(x) col_dend
      heatmap_color_scale <- function(data){
        data <- round(data, 3)
        if(max(data) <= 0 & min(data) < 0){
          over_median <- min(data)/2
          if(max(data) < over_median){
            color <- grDevices::colorRampPalette(
              c("#157AB5","#92c5de"))(n=1000)
          }else{
            color_rank <- round(max(data)/(min(data))*1000)
            color_scale <- grDevices::colorRampPalette(
              c("#0571b0","#92c5de","white"))(n=1000)
            color <- color_scale[color_rank:1000]
          }
        }else if(min(data) >= 0 & max(data) > 0){
          over_median <- max(data)/2
          if(min(data) > over_median){
            color <- grDevices::colorRampPalette(
              c("#f4a582", "#ca0020"))(n=1000)
          }else{
            color_rank <- round(min(data)/(max(data))*1000)
            color_scale <- grDevices::colorRampPalette(
              c("white","#f4a582", "#ca0020"))(n=1000)
            color <- color_scale[color_rank:1000]
          }
        }
        return(color)
      }
      if(min(data) >= 0 || max(data) <= 0){
        stats::heatmap(data, Rowv=TRUE, Colv=TRUE,
                       dendrogram='both', trace="none",
                       col=heatmap_color_scale(data),
                       distfun=sample_col_dend,
                       hclustfun=sample_col_dend,
                       main=NULL,
                       margins=c(8,8),
                       lwid=c(1, 9),
                       scale='none')
      }else{
        stats::heatmap(
          data, Rowv=TRUE, Colv=TRUE,
          dendrogram='both', trace="none",
          col=grDevices::colorRampPalette(c("blue",
                                            "white",
                                            "red"))(n=300),
          distfun=sample_col_dend,
          hclustfun=sample_col_dend,
          main=NULL,
          margins=c(8,8),
          lwid=c(1, 9),
          scale='none')
      }
      heatmap <- grDevices::recordPlot()
      grDevices::dev.off()
      reorder_data <- data[rev(col_dend$order),col_dend$order]
    }
  }else{
    heatmap <- NULL
    reorder_data <- NULL
  }
  
  return(list(corr_coef=corr_coef,
              corr_p=corr_p,
              heatmap=heatmap,
              reorder_data=reorder_data))
}
