#' @title Clin_LR_heatmap
#' @description Linear regression is a statistical technique that uses several explanatory variables to predict the outcome of a continuous response variable, allowing researchers to estimate the associations between lipid levels and clinical features. This function prints the heatmap based on Linear Regression results.
#' @param exp_data A data frame of predictors, including features (molecule, lipid class etc.) and their expression of each sample. NAs are allowed. The first column must be 'feature' of gene/lipid names.
#' @param condition_table A data frame. The condition table encompasses sample names and clinical conditions (disease status, gene dependence score, etc.)
#' @param adjusted_table A data frame with additional variables will be added into the algorithm and used to adjust the confounding effect.
#' @param adjust_p_method A character string of correction method. One of \bold{"holm"}, \bold{"hochberg"}, \bold{"hommel"}, \bold{"bonferroni"}, \bold{"BH"}, \bold{"BY"}, \bold{"fdr"}, \bold{"none"}, can be abbreviated. (default: "BH")
#' @param sig_stat A character string indicating which pvalue is to be used for the statistically significant. One of \bold{"p.adj"} or \bold{"p"}. (default: "p.adj")
#' @param sig_pvalue Numeric. Significant level.(default: 1)
#' @param distfun A character string of the distance measure indicating which correlation coefficient (or covariance) is to be computed. Allowed methods include \bold{"pearson"}, \bold{"kendall"}, \bold{"spearman"}.(default: "spearman")
#' @param hclustfun A character string of the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of \bold{"ward.D"}, \bold{"ward.D2"}, \bold{"single"}, \bold{"complete"}, \bold{"average"}(= UPGMA), \bold{"mcquitty"} (= WPGMA), \bold{"median"} (= WPGMC) or \bold{"centroid"} (= UPGMC). (default: "centroid")
#' @param heatmap_col A character string indicating the value for clustering. Allow method are \bold{'beta_coef'} and \bold{'t_statistic'}. (default: statistics)
#' @return Return a list with 2 data frames, 1 plot, and 1 matrix.
#' \enumerate{
#' \item Clin_LR_table_all: data frame showing statistical results from the correlation analysis.
#' \item Clin_LR_table_sig: data frame only showing significant statistical results from correlation analysis.
#' \item Clin_LR_table_plot: heatmap plot display correlation result(\bold{beta_coef'} or \bold{'t_statistic'}).
#' \item Clin_LR_reorder_mat: matrix of heatmap
#' }
#' @export
#' @examples
#' \dontrun{
#' data("corr_exp_data")
#' data("corr_condition_table")
#' data("corr_adjusted_table")
#' data("corr_adjusted_table")
#' exp_data <- corr_exp_data
#' condition_table <- corr_condition_table
#' adjusted_table <- corr_adjusted_table
#' exp_transform <- data_process(exp_data, exclude_var_missing=T, missing_pct_limit=50,
#'                               replace_zero=T, zero2what='min', xmin=0.5, replace_NA=T,
#'                               NA2what='min', ymin=0.5, pct_transform=T, data_transform=T,
#'                               trans_type='log', centering=F,  scaling=F)
#' Clin_LR_heatmap(exp_transform, condition_table, adjusted_table, adjust_p_method = 'BH',
#'                 sig_stat = 'p.adj', sig_pvalue = 1, distfun='spearman',
#'                 hclustfun='centroid', heatmap_col='beta_coef')
#' }
Clin_LR_heatmap <- function(exp_data,
                       condition_table, adjusted_table,
                       adjust_p_method = 'BH',
                       sig_stat = 'p.adj', sig_pvalue = 0.05,
                       distfun='spearman', hclustfun='centroid',
                       heatmap_col='beta_coef'){

  if(ncol(exp_data)<=10){
    stop("At least 10 samples.")
  }
  if(sum(sapply(exp_data[,-1], class)%in%c("numeric","integer"))!=ncol(exp_data[,-1])){
    stop("First column type must be 'character',others must be 'numeric'")
  }
  if(nrow(exp_data)!=length(unique(exp_data[,1]))){
    stop("The lipids name (features) must be unique")
  }
  if(sum(is.na(exp_data[,-1]))>0){
    stop("Variables can not be NA")
  }
  if(class(condition_table[,1])!="character"){
    stop("The first column must be 'sample_name'.")
  }
  if(ncol(condition_table)<3){
    stop("At least 2 conditions.")
  }
  if(nrow(condition_table)<3){
    stop("At least 10 samples.")
  }
  if(class(condition_table[,1])!="character" | sum(sapply(condition_table[,-1], class)%in%c("numeric","integer"))!=ncol(condition_table[,-1])){
    stop("The columns 'sample_name' must be characters; the other columns must be numeric.")
  }

  if(sum(condition_table[,1]%in%colnames(exp_data))!=nrow(condition_table) | sum(condition_table[,1]%in%colnames(exp_data))!=ncol(exp_data[,-1])){
    stop("'sample_name' must same as the name of samples of 'Lipid expression data'.")
  }
  if(class(adjusted_table[,1])!="character"){
    stop("The first column must contain a list of lipids names (features).")
  }
  if(sum(adjusted_table[,1]%in%colnames(exp_data))!=nrow(adjusted_table) | sum(adjusted_table[,1]%in%colnames(exp_data))!=ncol(exp_data[,-1])){
    stop("The sample names of 'Adjusted table' must same as the sample names of 'Lipid expression data' table.")
  }
  if(!distfun %in% c('pearson','spearman','kendall',"euclidean","maximum","manhattan","canberra","binary","minkowski")){
    stop('distfun must be one of the strings "pearson", "kendall", "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski"')
  }
  if(!hclustfun %in% c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")){
    stop('hclustfun must be one of the strings "ward.D","ward.D2","single","complete","average","mcquitty","median","centroid".')
  }

    CHAR <- colnames(exp_data)[1]
    colnames(exp_data)[1] <- 'feature'
    exp_data <- exp_data %>% tidyr::gather(-feature, key='sample_name', value='value') %>%
      tidyr::spread(key='feature', value='value') %>%
      dplyr::arrange(sample_name)

    condition_table <- condition_table %>% dplyr::arrange(sample_name)
    if(!is.null(adjusted_table)){
      adjusted_table <- adjusted_table %>% dplyr::arrange(sample_name)
    }

    # if(data_transform==T){
    #   exp_data[-1] <- log10(exp_data[-1])
    #   #condition_table[-1] <- log10(condition_table[-1])
    # }
    exp_data[-1] <- scale(exp_data[-1]) %>% as.data.frame()

    Clin_LR_table_all <- list(length = (ncol(condition_table)-1))
    Clin_LR_table_sig <- list(length = (ncol(condition_table)-1))
    #---------linear regression-----------------------------

    for(a in 2:ncol(condition_table)){
      feature <- character((ncol(exp_data)-1))
      beta_coef <- numeric((ncol(exp_data)-1))
      t_statistic <- numeric((ncol(exp_data)-1))
      p_value <- numeric((ncol(exp_data)-1))
      R_squared <- numeric((ncol(exp_data)-1))
      Adj_r_squared <- numeric((ncol(exp_data)-1))

      for(b in 2:ncol(exp_data)){

        feature[b-1] <- colnames(exp_data)[b]


        if(!is.null(adjusted_table)){
          data_lm <- cbind(condition_table[a],exp_data[b],adjusted_table[-1])
        }else{
          data_lm <- cbind(condition_table[a],exp_data[b])
        }
        colnames(data_lm)[1] <- c('clin_var')

        lm_summary <- tryCatch(
          summary(stats::lm(clin_var~., data=data_lm)),
          error=function(e){NULL})

        if(!is.null(lm_summary)){
          coef <- lm_summary$coefficients[2,]
          beta_coef[b-1] <- coef[1]
          t_statistic[b-1] <- coef[3]
          p_value[b-1] <- coef[4]
          R_squared[b-1] <- lm_summary$"r.squared"
          Adj_r_squared[b-1] <- lm_summary$"adj.r.squared"
        }else{
          beta_coef[b-1] <- NA
          t_statistic[b-1] <- NA
          p_value[b-1] <- NA
          R_squared[b-1] <- NA
          Adj_r_squared[b-1] <- NA
        }
      }

      Clin_LR_table_all[[a-1]] <- data.frame(clin_factor=colnames(condition_table)[a],
                                             feature=feature, method='Linear Regression',
                                             beta_coef=beta_coef, t_statistic=t_statistic,
                                             p_value=p_value, p_adj=stats::p.adjust(p_value, method = adjust_p_method, n = length(p_value)))

    }

    Clin_LR_table_all <- Reduce(rbind, Clin_LR_table_all)

    Clin_LR_table_all <- Clin_LR_table_all %>%
      dplyr::mutate(sig_p=ifelse(p_value<sig_pvalue, 'yes','no'),
             sig_p_adj=ifelse(p_adj<sig_pvalue, 'yes','no'))

    if(sig_stat=='p'){
      Clin_LR_table_sig <- Clin_LR_table_all %>%
        dplyr::filter(sig_p=='yes')
    }else if(sig_stat=='p.adj'){
      Clin_LR_table_sig <- Clin_LR_table_all %>%
        dplyr::filter(sig_p_adj=='yes')
    }


    #---------heatmap-----------------------------
    if(length(unique(Clin_LR_table_sig[[2]]))>1){
      max_colcex <- max(stringr::str_length(Clin_LR_table_sig[[2]]))
      max_rowcex <- max(stringr::str_length(Clin_LR_table_sig[[1]]))

      if(max_colcex<4){
        max_colcex <- 4
      }
      if(max_rowcex<4){
        max_rowcex <- 4
      }
      if(heatmap_col=='beta_coef'){
        Cor.mat <- Clin_LR_table_sig %>%
          dplyr::select(clin_factor, feature, beta_coef) %>%
          tidyr::spread(feature, beta_coef) %>%
          tibble::column_to_rownames(var = 'clin_factor') %>%
          as.matrix()
      }else if(heatmap_col=='t_statistic'){
        Cor.mat <- Clin_LR_table_sig %>%
          dplyr::select(clin_factor, feature, t_statistic) %>%
          tidyr::spread(feature, t_statistic) %>%
          tibble::column_to_rownames(var = 'clin_factor') %>%
          as.matrix()
      }


      Cor.mat[is.na(Cor.mat)] <- 0

      if(sum(is.na(exp_data))==0){
        cb_grid <- iheatmapr::setup_colorbar_grid(y_length =0.6,x_start = 1,y_start = 0.4)
        if(distfun%in%c("pearson","kendall","spearman")){
          col_dend <- stats::hclust(stats::as.dist(1-stats::cor(Cor.mat, method=distfun)),method = hclustfun)
          row_dend <- stats::hclust(stats::as.dist(1-stats::cor(t(Cor.mat), method=distfun)),method = hclustfun)
        }else{
          col_dend <- stats::hclust(stats::dist(t(Cor.mat), method=distfun),method = hclustfun)
          row_dend <- stats::hclust(stats::dist(Cor.mat, method=distfun),method = hclustfun)
        }

        heatmap_color_scale <- function(data){
          data <- round(data,3)
          if(max(data)<=0 & min(data)<0){
            over_median <- min(data)/2
            if(max(data)<over_median){
              color <-  grDevices::colorRampPalette(c("#157AB5","#92c5de"))(n = 1000)
            }else{
              color_rank <- round(max(data)/(min(data))*1000)
              color_scale <- grDevices::colorRampPalette(c("#0571b0","#92c5de","white"))(n = 1000)
              color <- color_scale[color_rank:1000]
            }
          }else if(min(data)>=0 & max(data)>0){
            over_median <- max(data)/2
            if(min(data)>over_median){
              color <-  grDevices::colorRampPalette(c("#f4a582", "#ca0020"))(n = 1000)
            }else{
              color_rank <- round(min(data)/(max(data))*1000)
              color_scale <- grDevices::colorRampPalette(c("white","#f4a582", "#ca0020"))(n = 1000)
              color <- color_scale[color_rank:1000]
            }
          }
          return(color)
        }
        #if(min(Cor.mat)>=0){
        #  hm <- iheatmapr::iheatmap(Cor.mat,colors = heatmap_color_scale(Cor.mat),colorbar_grid = cb_grid) %>%
        #    iheatmapr::add_col_dendro(col_dend,side="top",reorder =T,size=0.1) %>%
        #    iheatmapr::add_row_dendro(row_dend,side="right",reorder =T,size=0.1)
        #}else if(max(Cor.mat)<=0){
        #  hm <- iheatmapr::iheatmap(Cor.mat,colors = heatmap_color_scale(Cor.mat),colorbar_grid = cb_grid) %>%
        #    iheatmapr::add_col_dendro(col_dend,side="top",reorder =T,size=0.1) %>%
        #    iheatmapr::add_row_dendro(row_dend,side="right",reorder =T,size=0.1)
        #}else{
        #  hm <- iheatmapr::iheatmap(Cor.mat,colorbar_grid = cb_grid) %>%
        #    iheatmapr::add_col_dendro(col_dend,side="top",reorder =T,size=0.1) %>%
        #    iheatmapr::add_row_dendro(row_dend,side="right",reorder =T,size=0.1)
        #}
        #if(max(nchar(colnames(Cor.mat)))<10){
        #  text_size <- 0.1
        #}else if(max(nchar(colnames(Cor.mat)))<20 & max(nchar(colnames(Cor.mat)))>=10){
        #  text_size <- 0.2
        #}else if(max(nchar(colnames(Cor.mat)))<30 & max(nchar(colnames(Cor.mat)))>=20){
        #  text_size <- 0.3
        #}else if(max(nchar(colnames(Cor.mat)))<40 & max(nchar(colnames(Cor.mat)))>=30){
        #  text_size <- 0.4
        #}else{
        #  text_size <- 0.5
        #}
        #if(nrow(Cor.mat)<50){
        #  hm <- hm %>% iheatmapr::add_row_labels(font=list(size=8))
        #}
        #if(ncol(Cor.mat)<50){
        #  hm <- hm %>% iheatmapr::add_col_labels(side="bottom",font=list(size=8),size= text_size)
        #}
        ax <- list(
          title = "",
          zeroline = FALSE,
          showline = FALSE,
          showgrid = FALSE,
          ticks=''
        )
        bx <- list(
          title = "",
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = FALSE,
          ticks=''
        )
        if(distfun %in% c('pearson','spearman','kendall')){
          dist_fun=function(x){
            x=t(x)
            cor.mat=stats::cor(x,method=distfun,use = 'complete.obs')
            cor.mat=(1-cor.mat)
            cor.dist=stats::as.dist(cor.mat)
            return(cor.dist)
          }
        }else{
          dist_fun=function(x) stats::dist(x, method=distfun)
        }

        hclust_fun = function(x) stats::hclust(x, method = hclustfun)
        hm = Cor.mat %>%
          heatmaply::heatmaply(scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low="#0571b0",mid="white", high="#ca0020",midpoint = 0),
                               distfun = dist_fun,
                               hclustfun = hclust_fun,
                               column_text_angle = 270, grid_color = "black",
                               margins = c(l=0.2, r=0.8, t=20, b=80))
        if(nrow(Cor.mat)>50 & ncol(Cor.mat)>50){
          hm <- hm %>%
            plotly::layout(xaxis = bx, yaxis = bx)
        }else if(nrow(Cor.mat)>50){
          hm <- hm %>%
            plotly::layout(xaxis = ax, yaxis = bx)
        }else if(ncol(Cor.mat)>50){
          hm <- hm %>%
            plotly::layout(xaxis = bx, yaxis = ax)
        }else{
          hm <- hm %>%
            plotly::layout(xaxis = ax, yaxis = ax)
        }
        reorder_Cor.mat <- Cor.mat[rev(row_dend$order),col_dend$order]
      }else{
        hm = NULL
        reorder_Cor.mat = NULL
      }
    }else{
      hm = NULL
      reorder_Cor.mat = NULL
    }

    colnames(Clin_LR_table_all)[2] <- CHAR
    colnames(Clin_LR_table_sig)[2] <- CHAR

  return(list(LR_table_all=Clin_LR_table_all,
              LR_table_sig=Clin_LR_table_sig,
              LR_table_plot=hm,
              LR_reorder_mat=reorder_Cor.mat))
}
