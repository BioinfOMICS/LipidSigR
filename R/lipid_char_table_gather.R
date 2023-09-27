#' @title lipid_char_table_gather
#' @description Select the user-selected characteristics, and make it become 
#' a long table.
#' @param exp_data_SE A SummarizedExperiment object contains information about 
#' various features, such as molecules, lipid class, etc., and their 
#' expression levels for each sample. The assay is a matrix representing the 
#' expression of lipid features across all samples. Missing values (NAs)  are 
#' allowed. The row data  corresponds to specific lipid features, such as class 
#' and total length. The first column's name must be "feature" (lipid species), 
#' and NAs are allowed for this data. The column data comprises sample names, 
#' sample labels, group names, and pair numbers that represent 'the pair' for 
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of the rowData of \bold{exp_data_SE},
#' such as total length.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' library(SummarizedExperiment)
#' data("profiling_data")
#' exp_data_SE <- profiling_data
#' lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' lipid_char_gather <- lipid_char_table_gather(exp_data_SE, 
#'     char_var = char_var[8])
lipid_char_table_gather <- function(exp_data_SE, char_var){
  
  lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))

  if(!char_var %in% 
     grep(pattern='FA_', x=colnames(lipid_char_table), value=TRUE)){
    stop("Function lipid_char_table_gather is use for transform 'FA_' column")
  }
  lipid_char <- lipid_char_table %>%
    dplyr::select(feature, tidyselect::all_of(char_var))
  colnames(lipid_char) <- c('feature', 'characteristic')
  lipid_char_ga <- lipid_char %>%
    tidyr::separate(col=characteristic, into=c('c1', 'c2'), sep=',',
                    remove=TRUE, convert=TRUE, fill='right') %>%
    tidyr::gather(C, Value, -1) %>%
    dplyr::filter(!is.na(Value)) %>%
    dplyr::select(-C)
  colnames(lipid_char_ga) <- c('feature', char_var)
  
  lipid_char_ga_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(lipid_char_ga=lipid_char_ga))
  
  return(lipid_char_table_ga = lipid_char_ga_SE)
} 