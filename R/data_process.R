#' @title data_process
#' @description This function conducts data processing according to users'
#'     options, including removing features with missing values, missing values
#'     imputation, percentage transformation, log10 transformation, etc.
#' @param exp_data_SE A SummarizedExperiment object contains information about
#' various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data  corresponds to specific lipid features, such as class
#' and total length. The first column's name must be "feature" (lipid species),
#' and NAs are allowed for this data. The column data comprises sample names,
#' sample labels, group names, and pair numbers that represent 'the pair' for
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param exclude_var_missing Logical. Remove Lipid with too many
#' missing values. (default: TRUE)
#' @param missing_pct_limit An integer indicating the missing values over a
#' certain percentage should be removed. (default: 50)
#' @param replace_zero Logical. Replace 0. (default: TRUE)
#' @param zero2what A character string indicating the value to replace 0.
#' Value include "mean", "median", "min". (default: min)
#' @param xmin A numeric value indicating the min value to replace 0.
#' @param replace_NA Logical. If remove_na = TRUE, all NA will be removed.
#' (default: TRUE)
#' @param NA2what A character string indicating the value to replace NA.
#' Value include "mean", "median", "min". (default: min)
#' @param ymin A numeric value indicating the min value to replace NA.
#' @param pct_transform Logical. If \bold{pct_transform = TRUE}, transform
#' lipid value into a percentage. (default:TRUE)
#' @param data_transform Logical. If \bold{data_transform = TRUE},
#' transform exp_data by log10.
#' @param trans_type A character string of transformation type.
#' @param centering A logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal to
#' the number of columns of x can be supplied. The value is passed to scale.
#' (default: FALSE)
#' @param scaling A logical value. If scaling = TRUE, each block is
#' standardized to zero means and unit variances. (default: FALSE).
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' data("DE_data")
#' exp_data_SE <- DE_data
#' data_process(exp_data_SE, exclude_var_missing=TRUE, missing_pct_limit=50,
#'     replace_zero=TRUE, zero2what='min', xmin=0.5, replace_NA=TRUE,
#'     NA2what='min', ymin=0.5, pct_transform=TRUE, data_transform=TRUE,
#'     trans_type='log', centering=FALSE, scaling=FALSE)
data_process <- function(exp_data_SE, exclude_var_missing=TRUE,
                         missing_pct_limit=50, replace_zero=TRUE,
                         zero2what='min', xmin=0.5, replace_NA=TRUE,
                         NA2what='min', ymin=0.5, pct_transform=TRUE,
                         data_transform=TRUE,trans_type='log',
                         centering=FALSE, scaling=FALSE){
  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c(
    colnames(SummarizedExperiment::rowData(exp_data_SE))[[1]],
      SummarizedExperiment::colData(exp_data_SE)[[1]] )
  lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))

  exp_data2 <- exp_data[-1]
  min_exp <- min(unlist(exp_data2)[unlist(exp_data2)>0], na.rm=TRUE)
  if(replace_zero == TRUE){
    if(is.numeric(zero2what)){
      exp_data2[exp_data2 == 0] <- zero2what
    }else if(zero2what=='min'){
      exp_data2[exp_data2 == 0] <- min_exp*xmin
    }else if(zero2what=='NA'){
      exp_data2[exp_data2 == 0] <- NA
    }
    exp_data <- cbind(exp_data[1], exp_data2)
  }
  exp_data2 <- exp_data[-1]
  if(exclude_var_missing == TRUE){
    missing_pct <- apply(exp_data2, 1, function(x){sum(is.na(x)) / length(x)})
    maintain_var <- missing_pct*100 < missing_pct_limit
    exp_data <- exp_data[maintain_var,]
  }
  if(nrow(exp_data) == 0){
    return(NULL)
  }
  exp_data2 <- exp_data[-1]
  if(replace_NA == TRUE){
    if(NA2what == 'min'){
      for(a in seq_len(nrow(exp_data2))){
        exp_data2[a,][is.na(exp_data2[a,])] <- ymin*min_exp
      }
    }else if(NA2what == 'mean'){
      for(a in seq_len(nrow(exp_data2))){
        exp_data2[a,][is.na(exp_data2[a,])] <-
          mean(unlist(exp_data2[a,]), na.rm=TRUE)
      }
    }else if(NA2what == 'median'){
      for(a in seq_len(nrow(exp_data2))){
        exp_data2[a,][is.na(exp_data2[a,])] <-
          stats::median(unlist(exp_data2[a,]), na.rm=TRUE)
      }
    }else if(is.numeric(NA2what)){
      exp_data2[is.na(exp_data2)] <- NA2what
    }
    exp_data <- cbind(exp_data[1], exp_data2)
  }
  if(pct_transform == TRUE){
    exp_data[-1] <- purrr::map2(exp_data[-1],
      colSums(exp_data[-1], na.rm=TRUE), ~.x/.y*100) %>%
      as.data.frame()
  }
  if(data_transform == TRUE){
    if(trans_type == 'log'){exp_data[-1] <- log10(exp_data[-1])}
  }
  if(centering == TRUE){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale(scale=FALSE) %>% t() %>%
      as.data.frame()
  }
  if(scaling == TRUE){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale() %>% t() %>%
      as.data.frame()
  }

  lipid_char_table_trans <- as.data.frame(lipid_char_table[
    (lipid_char_table[[1]] %in% exp_data[[1]]), ])
  colnames(lipid_char_table_trans) <- colnames(lipid_char_table)

  transform_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(exp_data=as.matrix(exp_data[-1])),
    rowData=lipid_char_table_trans,
    colData=SummarizedExperiment::colData(exp_data_SE))

  return(transform_SE)
}
