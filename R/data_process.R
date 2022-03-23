#' @title data_process
#' @description This function conducts data processing according to users' options, including removing features with missing values, missing values imputation, percentage transformation, log10 transformation, etc.
#' @param exp_data A data frame of predictors, including features (molecule, lipid class, etc.) and their expression of each sample. NAs are not allowed. The name of the first column must be "feature" (lipid species).
#' @param exclude_var_missing Logical. Remove Lipid with too many missing values. (default: T)
#' @param missing_pct_limit An integer indicating the missing values over a certain percentage should be removed. (default: 50)
#' @param replace_zero Logical. Replace 0. (default: T)
#' @param zero2what A character string indicating the value to replace 0. Value include "mean", "median", "min". (default: min)
#' @param xmin A numeric value indicating the min value to replace 0.
#' @param replace_NA Logical. If remove_na = TRUE, all NA will be removed. (default: T)
#' @param NA2what A character string indicating the value to replace NA. Value include "mean", "median", "min". (default: min)
#' @param ymin A numeric value indicating the min value to replace NA.
#' @param pct_transform Logical. If \bold{pct_transform = TRUE}, transform lipid value into a percentage. (default:T)
#' @param data_transform Logical. If \bold{data_transform = TRUE}, transform exp_data by log10.
#' @param trans_type A character string of transformation type.
#' @param centering A logical value indicating whether the variables should be shifted to be zero centered. AAlternately, a vector of length equal to the number of columns of x can be supplied. The value is passed to scale. (default: F)
#' @param scaling A logical value. If scaling = TRUE, each block is standardized to zero means and unit variances. (default: F).
#' @return Return a list of 1 data frame.
#' @export
#' @examples
#' \donttest{
#' data("DE_exp_data")
#' exp_data <- DE_exp_data
#' data_process(exp_data, exclude_var_missing=T, missing_pct_limit=50, replace_zero=T,
#'              zero2what='min', xmin=0.5, replace_NA=T, NA2what='min', ymin=0.5,
#'              pct_transform=T, data_transform=T, trans_type='log', centering=F, scaling=F)
#' }
data_process <- function(exp_data, exclude_var_missing=T,
                         missing_pct_limit=50,
                         replace_zero=T, zero2what='min', xmin=0.5,
                         replace_NA=T, NA2what='min', ymin=0.5,
                         pct_transform=T,
                         data_transform=T,trans_type='log',
                         centering=F,
                         scaling=F){

  if(ncol(exp_data)==2){
    if(sum(class(exp_data[,-1])%in%c("numeric","integer"))!=1){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(sum(sapply(exp_data[,-1], class)%in%c("numeric","integer"))!=ncol(exp_data[,-1])){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }
  if(nrow(exp_data)!=length(unique(exp_data[,1]))){
    stop("The lipids name (features) must be unique")
  }
  if(ncol(exp_data)<3){
    stop("At least 2 samples.")
  }
  if(sum(exp_data[,-1]<0,na.rm = T)>0){
    stop("Variable must greater than zero")
  }
  if(sum(!is.na(exp_data[,-1]))==0 | sum(!is.null(exp_data[,-1]))==0){
    stop("Variables can not be all NULL/NA")
  }
  exp_data2 <- exp_data[-1]
  min_exp <- min(unlist(exp_data2)[unlist(exp_data2)>0],na.rm = T)
  #replace 0 with: min, specfic num
  if(replace_zero==T){
    if(is.numeric(zero2what)){
      exp_data2[exp_data2==0] <- zero2what
    }else if(zero2what=='min'){
      exp_data2[exp_data2==0] <- min_exp*xmin
    }else if(zero2what=='NA'){
      exp_data2[exp_data2==0] <- NA
    }
    exp_data <- cbind(exp_data[1],exp_data2)
  }



  exp_data2 <- exp_data[-1]
  #exclude_var_missing
  if(exclude_var_missing==T){
    missing_pct <- apply(exp_data2, 1, function(x){ sum(is.na(x)) / length(x) })
    maintain_var <- missing_pct*100 < missing_pct_limit
    exp_data <- exp_data[maintain_var,]
  }
  if(nrow(exp_data)==0){
    #stop('no species remains')
    return(NULL)
  }


  exp_data2 <- exp_data[-1]


  #replace NA with: min, mean, median, specfic num
  if(replace_NA==T){
    if(NA2what=='min'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- ymin*min_exp
      }
    }else if(NA2what=='mean'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- mean(unlist(exp_data2[a,]), na.rm = T)
      }
    }else if(NA2what=='median'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- stats::median(unlist(exp_data2[a,]), na.rm = T)
      }
    }else if(is.numeric(NA2what)){
      exp_data2[is.na(exp_data2)] <- NA2what
    }
    exp_data <- cbind(exp_data[1],exp_data2)
  }

  if(pct_transform==T){
    exp_data[-1] <- purrr::map2(exp_data[-1], colSums(exp_data[-1], na.rm = T), ~.x/.y*100) %>%
      as.data.frame()
  }

  #data_transform
  if(data_transform==T){
    if(trans_type=='log'){
      exp_data[-1] <- log10(exp_data[-1])
    }
  }

  #centering
  if(centering==T){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale(scale = F) %>% t() %>% as.data.frame()
  }

  #scaling
  if(scaling==T){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale() %>% t() %>% as.data.frame()
  }

  return(exp_data)
} #function: data_process()
