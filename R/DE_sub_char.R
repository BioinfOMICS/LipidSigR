#' @title DE_sub_char
#' @description This function conduct lipid Characteristics Analysis.
#' Lipid species are categorized and summarized into a new lipid expression
#' table according to two selected lipid characteristics, then conducted
#' differential expressed analysis.
#' @param exp_data_SE A SummarizedExperiment object contains information about
#' various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data  corresponds to specific lipid features, such as class
#' and total length. The first column's name must be "feature" (lipid species),
#' and NAs are allowed for this data. The column data comprises sample names,
#' sample labels, group names, and pair numbers that represent 'the pair' for
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param data_transform Logical. If data_transform = TRUE, transform exp_data
#' by log10.
#' @param split_var A character string of the second lipid characteristic
#' selected by users from the column name of the rowData of \bold{exp_data_SE}
#'  for the subgroup analysis, such as class. \emph{NOTE: This parameter will
#'  be used to split data before entering the main lipid characteristic '
#'  char_var'.}
#' @param char_var A character string of the first lipid characteristic
#' selected by users from tthe column name of the rowData of
#' \bold{exp_data_SE}, such as total length.
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
#' @return Return a list with 4 SummarizedExperiment objects.
#' \enumerate{
#' \item result_table1: scaled expression based on two user-selected lipid
#' characters' grouping
#' \item result_table2: two-way anova result with a post_hoc_test result.
#' \item result_table3: scaled expression based on the char_var's weighted mean.
#' \item result_table4: t.test that based on result_table3.
#' }
#' @export
#' @examples
#' data("DE_data")
#' exp_data_SE <- DE_data
#' lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' DE_sub_char(exp_data_SE, data_transform=TRUE, split_var = char_var[2],
#' char_var = char_var[4], paired = FALSE, sig_pvalue=0.05, sig_FC=2,
#' exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
#' zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
#' pct_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
DE_sub_char <- function(exp_data_SE, data_transform=TRUE, split_var, char_var,
      paired=FALSE, sig_pvalue=0.05, sig_FC=0,
      exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
      zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
      pct_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE){

  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature",
                          SummarizedExperiment::colData(exp_data_SE)[[1]])
  group_info <- as.data.frame(SummarizedExperiment::colData(exp_data_SE))
  rownames(group_info) <- NULL

  lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))
  if(ncol(
    dplyr::select(lipid_char_table,tidyselect::starts_with("FA_"))) == 0){
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
    if(methods::is(FA_lipid_char_table$lipid.category.value, 'character')){
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
      dplyr::filter(eval(parse(text=split_var)) == split_var_member[var_member])

    split_SE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(exp_data=split_exp_data[-1]),
      rowData=split_lipid_char_table,
      colData=group_info)

    spec2char <- Species2Char(split_SE, char_var)

    exp_data_trans <- data_process(spec2char, exclude_var_missing,
        missing_pct_limit, replace_zero, zero2what, xmin, replace_NA, NA2what,
        ymin, pct_transform, data_transform=FALSE, trans_type,
        centering, scaling)

    result_SE <- DE_char(exp_data_trans, data_transform=data_transform,
        paired=paired, sig_pvalue=sig_pvalue, sig_FC=sig_FC)

    result_1 <- SummarizedExperiment::assay(result_SE[[1]])
    result_2 <- SummarizedExperiment::assay(result_SE[[2]])
    result_3 <- SummarizedExperiment::assay(result_SE[[3]])
    result_4 <- SummarizedExperiment::assay(result_SE[[4]])

    result_1[split_var] <- split_var_member[var_member]
    result_2[split_var] <- split_var_member[var_member]
    result_3[split_var] <- split_var_member[var_member]
    result_4[split_var] <- split_var_member[var_member]
    result_table1[[var_member]] <- result_1 %>%
      dplyr::select(split_var, dplyr::everything())
    result_table2[[var_member]] <- result_2 %>%
      dplyr::select(split_var, dplyr::everything())
    result_table3[[var_member]] <- result_3 %>%
      dplyr::select(split_var, dplyr::everything())
    result_table4[[var_member]] <- result_4 %>%
      dplyr::select(split_var, dplyr::everything())
  }
  result_table1 <- Reduce(rbind, result_table1)
  result_table2 <- Reduce(rbind, result_table2)
  result_table3 <- Reduce(rbind, result_table3)
  result_table4 <- Reduce(rbind, result_table4)

  result_table1_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(result_table1=as.data.frame(result_table1)))
  result_table2_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(result_table2=as.data.frame(result_table2)))
  result_table3_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(result_table3=result_table3))
  result_table4_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(result_table4=as.data.frame(result_table4)))

  return(list(result_table1=result_table1_SE, result_table2=result_table2_SE,
              result_table3=result_table3_SE, result_table4=result_table4_SE))
}
