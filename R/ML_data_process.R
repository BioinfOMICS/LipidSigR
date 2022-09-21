#' @title ML_data_process
#' @description This function provides data preprocessing and output processed data for machine learning.
#' @param exp_data A data frame of predictors, including features (molecule, lipid class, etc.) and their expression of each sample. NAs are not allowed. The name of the first column must be "feature" (lipid species).
#' @param group_info A data frame of one clinical term (or one set of clinical terms) from demographic data.  The first column is the sample name and the second column is the group number assigned to each sample. The group number can only be 0 or 1. NAs are not allowed.
#' @param lipid_char_table A data frame of lipid characteristics such as name(feature) of lipid, class of lipid, the total length of lipid, and Fatty acid (FA_) related characteristics. NAs are allowed. The name of the first column must be "feature" (lipid species).
#' @param char_var A character string of the lipid characteristic selected by users from the column name of \bold{lipid_char_table}, such as 'class', 'total length'.
#' @param exclude_var_missing Logical. Remove Lipid with too many missing values. (default: TRUE)
#' @param missing_pct_limit An integer indicating the missing values over a certain percentage should be removed.
#' @param replace_zero Logical. Replace 0. (default: TRUE)
#' @param zero2what A character string indicating the value to replace 0. Value include 'mean', 'median', 'min'. (default: min)
#' @param xmin A numeric value indicating the min value to replace 0.
#' @param replace_NA Logical. If remove_na = TRUE, all NA will be removed. (default: TRUE)
#' @param NA2what A character string indicating the value to replace NA. Value include 'mean', 'median', 'min'. (default: min)
#' @param ymin A numeric value indicating the min value to replace NA.
#' @param pct_transform Logical. Transform lipid value into a percentage. (default:TRUE)
#' @param data_transform Logical. Data transformation(log). (default:TRUE)
#' @param trans_type A character string of data transformation type. (default:'log')
#' @param centering A logical value indicating whether the variables should be shifted to be zero centered. Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale. (default: FALSE)
#' @param scaling A logical value. If scaling = TRUE, each block is standardized to zero means and unit variances. (default: FALSE).
#' @return Return a list of 2 data frames.
#' \enumerate{
#' \item [[1]]: the processed data can be used as input of machine learning analysis. The value of this data frame is adding the lipid (feature)'s original expression to the expression of user-selected lipid characteristics.
#' \item [[2]]: The transpose data frame of \bold{[[1]]}, which is also processed data and can be used as the input of machine learning analysis. The value of this data frame is adding the lipid (feature)'s original expression to the expression of user-selected lipid characteristics.
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
#' ML_data_process(exp_data, group_info = condition_table, lipid_char_table,
#'                 char_var[1], exclude_var_missing=TRUE,
#'                 missing_pct_limit=50, replace_zero=TRUE,
#'                 zero2what='min', xmin=0.5,replace_NA=TRUE, NA2what='min',
#'                 ymin=0.5, pct_transform=TRUE, data_transform=TRUE,
#'                 trans_type='log', centering=FALSE, scaling=FALSE)
ML_data_process <- function(exp_data, group_info, lipid_char_table, char_var=NULL,
                              exclude_var_missing=TRUE,
                              missing_pct_limit=50,
                              replace_zero=TRUE, zero2what='min', xmin=0.5,
                              replace_NA=TRUE, NA2what='min', ymin=0.5,
                              pct_transform=TRUE,
                              data_transform=TRUE, trans_type='log',
                              centering=FALSE,
                              scaling=FALSE){
  
  exp_data <- as.data.frame(exp_data)
  group_info <- as.data.frame(group_info)
  lipid_char_table <- as.data.frame(lipid_char_table)
  if(!is(exp_data[,1], 'character')){
    stop("exp_data first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data)==2){
    if(!is(exp_data[,1], 'character') | sum(class(exp_data[,-1])%in%c("numeric","integer"))!=1){
      stop("exp_data first column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(!is(exp_data[,1], 'character') | sum(sapply(exp_data[,-1], class)%in%c("numeric","integer"))!=ncol(exp_data[,-1])){
      stop("exp_data first column type must be 'character',others must be 'numeric'")
    }
  }
  if(tibble::is.tibble(exp_data)){
    if(nrow(exp_data)!=nrow(unique(exp_data[,1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(exp_data)!=length(unique(exp_data[,1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if(ncol(exp_data)<61){
    stop("exp_data at least 60 samples.")
  }
  if(nrow(exp_data)<10){
    stop("exp_data number of lipids names (features) must be more than 10.")
  }
  if(sum(!is.na(exp_data[,-1]))==0 | sum(!is.null(exp_data[,-1]))==0){
    stop("exp_data variables can not be all NULL/NA")
  }
  if(ncol(group_info)==2){
    if(is(group_info[,1], 'character') & sum(group_info$group%in%c(0,1))==nrow(group_info)){
      if(sum(table(group_info[,2])>=30)!=2){
        stop("group_info each group must have more than 30 samples.")
      }
    }else{
      stop("group_info columns 'sample_name' must be characters; the column 'group' must be numeric (only be 0 or 1).")
    }
  }else{
    stop("group_info column must contain 'sample_name' (1st column), 'group' (2nd column)")
  }
  if(sum(group_info[,1]%in%colnames(exp_data))!=nrow(group_info) | sum(group_info[,1]%in%colnames(exp_data))!=ncol(exp_data[,-1])){
    stop("group_info column 'sample_name' must same as the name of samples of exp_data")
  }
  if(!is(lipid_char_table[,1], 'character')){
    stop("lipid_char_table first column must contain a list of lipids names (features).")
  }
  if(nrow(lipid_char_table)!=length(unique(lipid_char_table[,1]))){
    stop("lipid_char_table lipids names (features) must be unique.")
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
  if(!is.null(char_var)){
    if(!char_var %in% colnames(lipid_char_table)){
      stop("char_var must be included in the lipid_char_table.")
    }
  }
  if((!is.null(lipid_char_table))&&(!is.null(char_var))){
    char_num <- length(char_var)
    char_data <- list()
    for(a in seq_len(char_num)){
      char_data[[a]] <- Species2Char(exp_data, lipid_char_table, char_var[a]) %>% dplyr::mutate(type=char_var[a])
      char_data[[a]][[1]] <- stringr::str_c(char_var[a],'_',char_data[[a]][[1]])
      colnames(char_data[[a]])[1] <- 'feature'

    }
    exp_data2 <- exp_data %>% dplyr::mutate(type='species')
    data_raw <- rbind(Reduce(rbind, char_data), exp_data2)
    data <- data_raw[-ncol(data_raw)]
    rownames(data) <- data[[1]]
    data <- data[-1] %>% t() %>% as.data.frame() %>%
      dplyr::mutate(sample_name=colnames(exp_data)[-1])
  }else{
    data_raw <- exp_data
    data_raw <- data_raw %>% dplyr::mutate(type='species')
    rownames(exp_data) <- exp_data[[1]]
    data <- exp_data[-1] %>% t() %>% as.data.frame() %>%
      dplyr::mutate(sample_name=colnames(exp_data)[-1])
  }

  data <- data %>% dplyr::left_join(group_info, by='sample_name') %>%
    dplyr::select(-sample_name) %>% dplyr::select(group, dplyr::everything())

  ML_data <- list(data_raw, data)

  combo <- unique(ML_data[[1]]$type)
  char_data <- list()
  for(var_num in seq_len(length(combo))){
    trans_data <- ML_data[[1]] %>% dplyr::filter(type==combo[var_num]) %>% dplyr::select(-type)
    char_data[[var_num]] <- data_process(trans_data, exclude_var_missing=exclude_var_missing,
                                         missing_pct_limit=missing_pct_limit,
                                         replace_zero=FALSE, zero2what, xmin,
                                         replace_NA=replace_NA, NA2what=NA2what, ymin=ymin,
                                         pct_transform=pct_transform,
                                         data_transform=FALSE, trans_type=FALSE,
                                         centering=FALSE, scaling=FALSE)
  }
  ML_data[[1]] <- Reduce(rbind, char_data)
  ML_data[[1]] <- data_process(ML_data[[1]], exclude_var_missing=FALSE,
                               missing_pct_limit=missing_pct_limit,
                               replace_zero=FALSE, zero2what, xmin,
                               replace_NA=FALSE, NA2what=NA2what, ymin=ymin,
                               pct_transform=FALSE,
                               data_transform, trans_type,
                               centering, scaling)
  data <- ML_data[[1]]
  rownames(data) <- data$feature
  ML_data[[2]] <- data[-1] %>% t() %>% as.data.frame() %>% dplyr::mutate(group=ML_data[[2]]$group) %>%
    dplyr::select(group, dplyr::everything())
  return(ML_data)
}
