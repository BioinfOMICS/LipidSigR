#' @title DE_sub_char_2
#' @description In the \code{\link{DE_sub_char_2}}, lipid Characteristics
#' Analysis will be conducted. Lipid species are categorized and summarized
#' into a new lipid expression table according to two selected
#' lipid characteristics, then conducted differential expressed analysis. \\
#' Note:
#' \enumerate{
#' \item The subgroup analysis is based on the first user-selected
#' characteristics. So, use \code{\link{DE_char_2}} before using
#' \code{\link{DE_sub_char_2}},
#' \item two selected characteristics should be both continuous data and one
#' categorical data with one continuous data. The cut-offs of differentially
#' expressed lipids are inputted by users.
#' }
#' @param exp_data A data frame includes the expression of lipid features in
#' each sample. NAs are allowed. First column should be gene/lipid name and
#' first column name must be 'feature'.
#' @param data_transform Logical. If data_transform = TRUE, transform exp_data
#' by log10.
#' @param lipid_char_table A data frame with lipid features, such as class,
#' total length. NAs are allowed. The name of the first column must be
#' "feature".
#' @param split_var A character string of the second lipid characteristic
#' selected by users from \bold{lipid_char_table} for the subgroup analysis,
#' such as class. \emph{NOTE: This parameter will be used to split data before
#' entering the main lipid characteristic 'char_var'.}
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of \bold{lipid_char_table},
#' such as total length.
#' @param group_info A data frame comprises the name of the sample,
#' the label of the sample, the group name of the sample, and the pair number
#' represents 'the pair' for the t-test/Wilcoxon test. NAs are allowed.
#' @param paired Logical. If paired = TRUE, data are paired samples.
#' (default: FALSE)
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @param sig_FC Numeric. Significance of the fold-change. (default: 0)
#' @param exclude_var_missing Logical. Remove Lipid with too many missing
#' values. (default: TRUE)
#' @param missing_pct_limit An integer indicating the missing values over a
#' certain percentage should be removed. (default: 50)
#' @param replace_zero Logical. Replace 0. (default: TRUE)
#' @param zero2what A character string indicating the value to replace 0. Value
#' include "mean", "median", "min". (default: min)
#' @param xmin A numeric value indicating the min value to replace 0.
#' @param replace_NA Logical. If remove_na = TRUE, all NA will be removed.
#' (default: TRUE)
#' @param NA2what A character string indicating the value to replace NA. Value
#' include "mean", "median", "min". (default: min)
#' @param ymin A numeric value indicating the min value to replace NA.
#' @param pct_transform Logical. If \bold{pct_transform = TRUE}, transform
#' lipid value into a percentage. (default:TRUE)
#' @param trans_type Logical. If \bold{data_transform = TRUE}, transform
#' exp_data by log10.
#' @param centering A logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal to the
#' number of columns of x can be supplied. The value is passed to scale.
#' (default: FALSE)
#' @param scaling A logical value. If scaling = TRUE, each block is
#' standardized to zero means and unit variances. (default: FALSE).
#' @return Return a list with 4 data frames.
#' \enumerate{
#' \item result_table1: a data frame with scaled expression based on two
#' user-selected lipid characters' grouping
#' \item result_table2: a data frame with two-way anova result with a
#' post_hoc_test result
#' \item result_table3: a data frame of scaled expression based on the
#' char_var's weighted mean
#' \item result_table4: a data frame of t.test that based on result_table3
#' }
#' @export
#' @examples
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' group_info <- DE_group_info[c(1:5, 11:15), ]
#' exp_data <- DE_exp_data[1:100, c("feature", group_info$sample_name)]
#' lipid_char_table <- DE_lipid_char_table[1:100, ]
#' char_var <- colnames(lipid_char_table)[-1]
#' DE_sub_char_2(exp_data, data_transform=TRUE,
#'               lipid_char_table = lipid_char_table,
#'               split_var = char_var[2], char_var = char_var[6],
#'               group_info = group_info,
#'               paired = FALSE, sig_pvalue=0.05, sig_FC=2,
#'               exclude_var_missing=TRUE, missing_pct_limit=50,
#'               replace_zero=TRUE, zero2what='min', xmin=0.5,
#'               replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
#'               trans_type='log', centering=FALSE, scaling=FALSE)
DE_sub_char_2 <- function(exp_data, data_transform=TRUE, lipid_char_table,
                          split_var, char_var,
                          group_info, paired=FALSE, sig_pvalue=0.05, sig_FC=0,
                          exclude_var_missing=TRUE, missing_pct_limit=50,
                          replace_zero=TRUE, zero2what='min', xmin=0.5,
                          replace_NA=TRUE, NA2what='min', ymin=0.5,
                          pct_transform=TRUE,
                          trans_type='log',
                          centering=FALSE, scaling=FALSE){

  exp_data <- as.data.frame(exp_data)
  if(!is(exp_data[,1], 'character')){
    stop("exp_data first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data)==2){
    if(!is(exp_data[, 1], 'character') |
       sum(class(exp_data[, -1]) %in% c("numeric", "integer")) != 1){
      stop("exp_data first column type must be 'character',
           others must be 'numeric'")
    }
  }else{
    if(!is(exp_data[,1], 'character') |
       sum(vapply(exp_data[, -1], class, character(1)) %in%
           c("numeric","integer")) != ncol(exp_data[, -1])){
      stop("exp_data first column type must be 'character',others must be
           'numeric'")
    }
  }
  if(tibble::is_tibble(exp_data)){
    if(nrow(exp_data) != nrow(unique(exp_data[, ]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(exp_data) != length(unique(exp_data[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if(ncol(exp_data) < 3){
    stop("exp_data at least 2 samples.")
  }else if(ncol(exp_data) == 3){
    warning("exp_data only 2 samples will not show p-value,dotchart will
            color by log2FC")
  }
  if(nrow(exp_data)<5){
    stop("exp_data number of lipids names (features) must be more than 5.")
  }
  if(sum(!is.na(exp_data[,-1]))==0 | sum(!is.null(exp_data[,-1]))==0){
    stop("exp_data variables can not be all NULL/NA")
  }
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
       sum(table(group_info[, 4]) != 2) != 0 & sum(is.na(group_info[, 4])) != 0)
      {
      stop("group_info each pair must have a specific number, staring from
           1 to N. Cannot have NA, blank, or skip numbers.")
    }
    if(sum(group_info[, 1] %in% colnames(exp_data)) != nrow(group_info) |
       sum(group_info[, 1] %in% colnames(exp_data)) != ncol(exp_data[, -1])){
      stop("group_info 'sample_name' must same as the name of
           samples of exp_data")
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
  }else if(ncol(group_info) == 3){
    if("pair" %in% colnames(group_info)){
      stop("group_info column must arrange in order of sample_name,
           label_name, group, pair(optional).")
    }
    if(sum(vapply(group_info,class,character(1)) != "character") != 0){
      stop("group_info first 3 columns must be characters.")
    }
    if(sum(group_info[, 1] %in% colnames(exp_data)) != nrow(group_info) |
       sum(group_info[, 1] %in% colnames(exp_data)) != ncol(exp_data[, -1])){
      stop("group_info 'sample_name' must same as the name of samples
           of exp_data")
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
  }
  if(nrow(lipid_char_table) == nrow(exp_data)){
    if(sum(lipid_char_table[, 1] %in% exp_data[, 1]) != nrow(lipid_char_table)){
      stop("The lipids names (features) of lipid_char_table table must
           same as exp_data.")
    }
  }else{
    stop("The row number of lipid_char_table table must same as exp_data.")
  }
  if(!is(lipid_char_table[, 1], 'character')){
    stop("lipid_char_table first column must contain a list of
         lipids names (features).")
  }
  if(nrow(lipid_char_table) != length(unique(lipid_char_table[, 1]))){
    stop("lipid_char_table lipids names (features) must be unique.")
  }
  if("class" %in% colnames(lipid_char_table)){
    if(!is(lipid_char_table[, 'class'], 'character')){
      stop("lipid_char_table content of column 'class' must be characters")
    }
  }
  if("totallength" %in% colnames(lipid_char_table)){
    if(!class(lipid_char_table[, 'totallength']) %in% c("integer", "numeric")){
      stop("lipid_char_table content of column 'totallength' must be numeric")
    }
  }
  if("totaldb" %in% colnames(lipid_char_table)){
    if(!class(lipid_char_table[, 'totaldb']) %in% c("integer", "numeric")){
      stop("lipid_char_table content of column 'totaldb' must be numeric")
    }
  }
  if("totaloh" %in% colnames(lipid_char_table)){
    if(!class(lipid_char_table[, 'totaloh']) %in% c("integer", "numeric")){
      stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
    }
  }

  if(ncol(dplyr::select(lipid_char_table,tidyselect::starts_with("FA_"))) == 0){
    warning("(OPTIONAL) lipid_char_table does not contain column names
            starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>%
      dplyr::select(feature, tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_", colnames(FA_lipid_char_table), value=TRUE)
    max_comma <- 0
    for(i in seq_len(length(FA_col))){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[,col], ','),
                         na.rm = TRUE)
      if(comma_count>0){
        FA_lipid_char_table <- FA_lipid_char_table %>%
          tidyr::separate(col, c(col,paste0(col, "_", seq_len(comma_count))),
                          ",", convert = TRUE)
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>%
      tidyr::gather(lipid.category, lipid.category.value,-feature)
    if(max_comma>0){
      for (i in seq_len(max_comma)) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <- FA_lipid_char_table[-intersect(
          grep(select_name, FA_lipid_char_table[, "lipid.category"]),
          which(is.na(FA_lipid_char_table$lipid.category.value))), ]
      }
    }
    if(is(FA_lipid_char_table$lipid.category.value, 'character')){
      stop("In the 'FA_' related analyses, the values are positive integer or
           zero and separated by comma. i.e., 10,12,11")
    }else if(sum(
      stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)) !=
      round(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value)))) != 0 |
      min(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value))) < 0){
      stop("In the 'FA_' related analyses, the values are positive integer or
           zero and separated by comma. i.e., 10,12,11")
    }
  }
  if(!char_var %in% colnames(lipid_char_table)){
    stop("char_var must be included in the lipid_char_table.")
  }
  split_var_member <- lipid_char_table %>%
    dplyr::filter(!is.na(eval(parse(text=char_var)))) %>%
    .[[split_var]] %>% unique() %>% stats::na.omit()

  exp_data <- exp_data %>%
    dplyr::left_join(lipid_char_table[c('feature', split_var)], by='feature')

  result_table1 <- list()
  result_table2 <- list()
  result_table3 <- list()
  result_table4 <- list()

  for(var_member in seq_len(length(split_var_member))){
    split_exp_data <- exp_data %>%
      dplyr::filter(eval(parse(text=split_var)) ==
                      split_var_member[var_member]) %>%
      dplyr::select(-split_var)
    split_lipid_char_table <- lipid_char_table %>%
      dplyr::filter(eval(parse(text=split_var)) ==
                      split_var_member[var_member])

    spec2char <- Species2Char(split_exp_data, split_lipid_char_table, char_var)

    if(nrow(spec2char)==0){
      methods::show('no')
      next()
    }

    exp_data_trans <- LipidSigR::data_process(spec2char, exclude_var_missing,
                                              missing_pct_limit,
                                              replace_zero, zero2what, xmin,
                                              replace_NA, NA2what, ymin,
                                              pct_transform,
                                              data_transform=FALSE, trans_type,
                                              centering, scaling)

    if(is.null(exp_data_trans)){
      methods::show('no')
      next()
    }

    result <- DE_char_2(exp_data_trans, data_transform=data_transform,
                        group_info=group_info, paired=paired,
                        sig_pvalue=sig_pvalue, sig_FC=sig_FC)

    result[[1]][split_var] <- split_var_member[var_member]
    result[[2]][split_var] <- split_var_member[var_member]
    result[[3]][split_var] <- split_var_member[var_member]
    result[[4]][split_var] <- split_var_member[var_member]

    result_table1[[var_member]] <- result[[1]] %>%
      dplyr::select(split_var, dplyr::everything())
    result_table2[[var_member]] <- result[[2]] %>%
      dplyr::select(split_var, dplyr::everything())
    result_table3[[var_member]] <- result[[3]] %>%
      dplyr::select(split_var, dplyr::everything())
    result_table4[[var_member]] <- result[[4]] %>%
      dplyr::select(split_var, dplyr::everything())


  }
  result_table1 <- Reduce(rbind, result_table1)
  result_table2 <- Reduce(rbind, result_table2)
  result_table3 <- Reduce(rbind, result_table3)
  result_table4 <- Reduce(rbind, result_table4)

  return(list(result_table1, result_table2, result_table3, result_table4))
}
