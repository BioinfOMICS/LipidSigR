#' @title lipid_char_table_gather
#' @description Select the user-selected characteristics, and make it become 
#' a long table.
#' @param lipid_char_table A data frame of lipid characteristics such as 
#' name(feature) of lipid, class of lipid, the total length of lipid, and 
#' Fatty acid (FA_) related characteristics. NAs are allowed. The name of the 
#' first column must be "feature" (lipid species).
#' @param char_var A character string of the lipid characteristic selected by 
#' users from the column name of \bold{lipid_char_table}, such as 'class'.
#' @return Return 1 long matrix of lipid species with user-selected 
#' lipid characteristics.
#' @export
#' @examples
#' data("profiling_lipid_char_table")
#' lipid_char_table <- profiling_lipid_char_table
#' char_var <- colnames(lipid_char_table)[-1]
#' lipid_char_table_gather(lipid_char_table, char_var = char_var[8])
lipid_char_table_gather <- function(lipid_char_table, char_var){
  
  lipid_char_table <- as.data.frame(lipid_char_table)
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
  
  return(lipid_char_table_ga = lipid_char_ga)
  
} #function
