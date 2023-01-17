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
#' @return Return 2 plots and a list of 6 matrices.
#' \enumerate{
#' \item sample_corr_coef: a matrix of correlation coefficients between samples
#' \item sample_corr_p: a matrix of correlation p-value between samples
#' \item reorder_sample_corr_coef: a matrix of reorder sample correlation.
#' \item lipids_corr_coef: a matrix of correlation coefficients
#' between features (lipids).
#' \item lipids_corr_p: a matrix of correlation p-value
#' between features (lipids).
#' \item reorder_lipids_corr_coef: reorder features (lipids)
#' correlation coefficients matrix
#' \item heatmap plot by samples
#' \item heatmap plot by lipid species
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
#'              hclustfun="average", plotly=TRUE)
corr_heatmap <- function(exp_data,corr_method="pearson",
                         distfun="maximum", hclustfun='average',
                         plotly=TRUE){
  exp_data <- as.data.frame(exp_data)
  #####################################
  #####################################
  ### Correlation heatmap by Sample ###
  #####################################
  #####################################
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
  ## Sample correlation coefficients ##
    sample_exp_data <- subset(exp_data, select=-feature)
    sample_corr_coef <- stats::cor(sample_exp_data,
                                   method=corr_method,
                                   use="pairwise.complete.obs")

    ## Sample correlation p-value ##
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
    sample_corr_p <- cor.test.p(x=sample_exp_data, method=corr_method)

    heatmap_color_scale <- function(data){
      data <- round(data,3)
      if(max(data) <= 0 & min(data) < 0){
        over_median <- min(data)/2
        if(max(data) < over_median){
          color <-  grDevices::colorRampPalette(c("#157AB5", "#92c5de"))(n=1000)
        }else{
          color_rank <- round(max(data)/(min(data))*1000)
          color_scale <- grDevices::colorRampPalette(c("#0571b0",
                                                       "#92c5de",
                                                       "white"))(n=1000)
          color <- color_scale[color_rank:1000]
        }
      }else if(min(data) >= 0 & max(data) > 0){
        over_median <- max(data)/2
        if(min(data)>over_median){
          color <-  grDevices::colorRampPalette(c("#f4a582", "#ca0020"))(n=1000)
        }else{
          color_rank <- round(min(data)/(max(data))*1000)
          color_scale <- grDevices::colorRampPalette(c("white",
                                                       "#f4a582",
                                                       "#ca0020"))(n=1000)
          color <- color_scale[color_rank:1000]
        }
      }
      return(color)
    }

    ## Sample Correlation heatmap ##
    cb_grid <- iheatmapr::setup_colorbar_grid(y_length=0.6,
                                              x_start=1,
                                              y_start=0.5)

    if(sum(is.na(sample_corr_coef)) == 0 &
       nrow(sample_corr_coef) >= 2 &
       ncol(sample_corr_coef) >= 2){
      if(distfun %in% c("pearson", "kendall", "spearman")){
        col_dend <- stats::hclust(
          stats::as.dist(1-stats::cor(sample_exp_data, method=distfun)),
          method=hclustfun)
      }else{
        col_dend <- stats::hclust(
          stats::dist(t(sample_exp_data), method=distfun), method=hclustfun)
      }
      sample_col_dend <- function(x) col_dend
      if(plotly == TRUE){
        if(min(sample_corr_coef) >= 0 || max(sample_corr_coef) <= 0){
          sample_hm <- iheatmapr::iheatmap(
            sample_corr_coef, colors=heatmap_color_scale(sample_corr_coef),
            colorbar_grid=cb_grid)
        }else{
          sample_hm <- iheatmapr::iheatmap(sample_corr_coef,
                                           colorbar_grid=cb_grid)
        }
        sample_hm <- sample_hm %>%
          iheatmapr::add_col_dendro(col_dend, side="top",
                                    reorder=TRUE, size=0.1) %>%
          iheatmapr::add_row_dendro(col_dend, side="right",
                                    reorder=TRUE, size=0.1)
        if(nrow(sample_corr_coef) < 50){
          sample_hm <- sample_hm %>%
            iheatmapr::add_row_labels(font=list(size=10))
        }
        if(ncol(sample_corr_coef) < 50){
          sample_hm <- sample_hm %>%
            iheatmapr::add_col_labels(side="bottom", font=list(size=10))
        }
        ## reorder Sample correlation matrix ##
        reorder_sample_corr_coef <- apply(sample_corr_coef[, col_dend$order],
                                          2, rev)
      }else{
        if(min(sample_corr_coef) >= 0 || max(sample_corr_coef) <= 0){
          stats::heatmap(sample_corr_coef, Rowv=TRUE, Colv=TRUE,
                         dendrogram='both', trace="none",
                         col=heatmap_color_scale(sample_corr_coef),
                         distfun=sample_col_dend,
                         hclustfun=sample_col_dend,
                         main=NULL,
                         margins=c(8,8),
                         lwid=c(1, 9),
                         scale='none')
        }else{
          stats::heatmap(
            sample_corr_coef, Rowv=TRUE, Colv=TRUE,
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
        sample_hm <- grDevices::recordPlot()
        grDevices::dev.off()
        reorder_sample_corr_coef <- sample_corr_coef[rev(col_dend$order),
                                                     col_dend$order]
      }
    }else{
      sample_hm <- NULL
      reorder_sample_corr_coef <- NULL
    }

    #####################################
    #####################################
    ### Correlation heatmap by lipids ###
    #####################################
    #####################################

    ## Lipids correlation coefficients ##
    lipids_exp_data <- exp_data %>%
      tidyr::gather(-feature, key='sample_name', value='value') %>%
      tidyr::spread(key='feature', value='value')
    lipids_exp_data <- subset(lipids_exp_data, select=-sample_name)
    lipids_corr_coef <- stats::cor(lipids_exp_data, method=corr_method,
                                   use="pairwise.complete.obs")

    ## Lipids correlation p-value ##
    lipids_corr_p <- cor.test.p(x=lipids_exp_data, method=corr_method)

    ## Lipids Correlation heatmap ##
    if(sum(is.na(lipids_exp_data)) == 0 &
       nrow(lipids_exp_data) >= 2 & ncol(lipids_exp_data) >= 2){
      if(distfun %in% c("pearson", "kendall", "spearman")){
        col_dend <- stats::hclust(
          stats::as.dist(1-stats::cor(lipids_exp_data, method=distfun)),
          method=hclustfun)
      }else{
        col_dend <- stats::hclust(
          stats::dist(t(lipids_exp_data), method=distfun), method=hclustfun)
      }
      lipids_col_dend <- function(x) col_dend
      if(plotly == TRUE){
        if(min(lipids_corr_coef) > 0 | max(lipids_corr_coef) < 0){
          lipids_hm <- iheatmapr::iheatmap(
            lipids_corr_coef, colors=heatmap_color_scale(lipids_corr_coef),
            colorbar_grid=cb_grid)
        }else{
          lipids_hm <- iheatmapr::iheatmap(
            lipids_corr_coef,colorbar_grid=cb_grid)
        }
        lipids_hm <- lipids_hm %>%
          iheatmapr::add_col_dendro(col_dend, side="top",
                                    reorder =TRUE, size=0.1) %>%
          iheatmapr::add_row_dendro(col_dend, side="right",
                                    reorder =TRUE, size=0.1)
        if(nrow(lipids_corr_coef) < 50){
          lipids_hm <- lipids_hm %>%
            iheatmapr::add_row_labels(font=list(size=10))
        }
        if(ncol(lipids_corr_coef) < 50){
          lipids_hm <- lipids_hm %>%
            iheatmapr::add_col_labels(side="bottom", font=list(size=10))
        }
        ## reorder Lipids correlation matrix ##
        reorder_lipids_corr_coef <- apply(lipids_corr_coef[,col_dend$order],
                                          2, rev)
      }else{
        if(min(lipids_corr_coef) > 0 | max(lipids_corr_coef) < 0){
          stats::heatmap(
            lipids_corr_coef, Rowv=TRUE, Colv=TRUE,
            dendrogram='both', trace="none",
            col=heatmap_color_scale(lipids_corr_coef),
            distfun=lipids_col_dend,
            hclustfun=lipids_col_dend,
            main=NULL,
            margins=c(8,8),
            lwid=c(1, 9),
            scale='none')
        }else{
          stats::heatmap(
            lipids_corr_coef, Rowv=TRUE, Colv=TRUE,
            dendrogram='both', trace="none",
            col=grDevices::colorRampPalette(c("blue",
                                              "white",
                                              "red"))(n=300),
            distfun=lipids_col_dend,
            hclustfun=lipids_col_dend,
            main=NULL,
            margins=c(8,8),
            lwid=c(1, 9),
            scale='none')
        }
        lipids_hm <- grDevices::recordPlot()
        grDevices::dev.off()
        reorder_lipids_corr_coef <- lipids_corr_coef[rev(col_dend$order),
                                                     col_dend$order]
      }
    }else{
      lipids_hm <- NULL
      reorder_lipids_corr_coef <- NULL
    }


  return(list(sample_corr_coef=sample_corr_coef,
              sample_corr_p=sample_corr_p,
              sample_hm=sample_hm,
              reorder_sample_corr_coef=reorder_sample_corr_coef,
              lipids_corr_coef=lipids_corr_coef,
              lipid_corr_p=lipids_corr_p,
              lipids_hm=lipids_hm,
              reorder_lipids_corr_coef=reorder_lipids_corr_coef))
}





