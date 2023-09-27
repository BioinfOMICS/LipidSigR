#' @title Species2Char
#' @description Summary the expression in \bold{exp_data_SE} from lipid species
#' to user-selected characteristics.
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
#' data("DE_data")
#' exp_data_SE <- DE_data
#' lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' Species2Char(exp_data_SE, char_var=char_var[1])
Species2Char <- function(exp_data_SE, char_var){

  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature",
                          SummarizedExperiment::colData(exp_data_SE)[[1]])

  lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))

  exp_data[-1][is.na(exp_data[-1])] <- 0
  if(stringr::str_detect(char_var, 'FA_')){
    max_chain_num <- stringr::str_split(lipid_char_table[[char_var]],
                                        pattern = ',') %>%
      purrr::map_dbl(length) %>% max()
    transform_table <- exp_data %>%
      dplyr::left_join(lipid_char_table[c('feature', char_var)],
                       by='feature') %>%
      tidyr::separate(eval(parse(text = char_var)),
                      letters[seq_len(max_chain_num)]) %>%
      tidyr::gather(letters[seq_len(max_chain_num)],
                    key='key', value='value') %>%
      dplyr::filter(!is.na(value))
    if(nrow(transform_table) == 0){
      transform_table <- data.frame()
    }else{
      transform_table <- transform_table %>%
        dplyr::select(-feature,-key) %>%
        stats::aggregate(. ~ value, ., sum)
      transform_table[[1]] <- as.character(transform_table[[1]])
      colnames(transform_table)[1] <- char_var
    }
  }else{
    transform_table <- exp_data %>%
      dplyr::left_join(lipid_char_table[c('feature',
                                          char_var)], by='feature') %>%
      dplyr::select(-1)
    transform_table <- transform_table[!is.na(transform_table[[char_var]]),]
    if(nrow(transform_table) == 0){
      transform_table <- data.frame()
    }else{
      transform_table <- transform_table %>%
        stats::aggregate(
          stats::as.formula(stringr::str_c('. ~ ',char_var)), ., sum)
    }
  }
  transform_table[-1][transform_table[-1] == 0] <- NA

  transform_mat <- as.matrix(transform_table[-1])
  colnames(transform_mat) <- NULL
  lipid_char_table_trans <- transform_table[1]

  transform_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(transform_table=transform_mat),
    rowData=lipid_char_table_trans,
    colData=SummarizedExperiment::colData(exp_data_SE))

  return(transform_SE)
}
