#' @title Species2Char
#' @description Summary the expression in \bold{exp_data} from lipid species
#' to user-selected characteristics.
#' @param exp_data A data frame of predictors, including features
#' (molecule, lipid class, etc.) and their expression of each sample. NAs are
#' not allowed. The name of the first column must be "feature" (lipid species).
#' @param lipid_char_table A data frame of lipid characteristics such as
#' name(feature) of lipid, class of lipid, the total length of lipid, and
#' Fatty acid (FA_) related characteristics. NAs are allowed. The name of the
#' first column must be "feature" (lipid species).
#' @param char_var A character string of the lipid characteristic selected by
#' users from the column name of \bold{lipid_char_table}, such as 'class'
#' @return Return \bold{transform_table}. An data frame of aggregated(sum)
#' expression data by user-selected characteristics (\bold{char_var}).
#' @export
#' @examples
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' char_var <- colnames(lipid_char_table)[-1]
#' Species2Char(exp_data, lipid_char_table, char_var=char_var[1])
Species2Char <- function(exp_data, lipid_char_table, char_var){

  exp_data <- as.data.frame(exp_data)
  lipid_char_table <- as.data.frame(lipid_char_table)
  if(nrow(lipid_char_table) == nrow(exp_data)){
    if(sum(lipid_char_table[,1] %in% exp_data[,1]) != nrow(lipid_char_table)){
      stop("The lipids names (features) of lipid_char_table
           table must same as exp_data.")
    }
  }else{
    stop("The row number of lipid_char_table table must same as exp_data.")
  }
  if(!is(lipid_char_table[,1], 'character')){
    stop("lipid_char_table first column must contain a list of
         lipids names (features).")
  }
  if(tibble::is_tibble(lipid_char_table)){
    if(nrow(lipid_char_table) != nrow(unique(lipid_char_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(lipid_char_table)!=length(unique(lipid_char_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
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
    if(!class(lipid_char_table[, 'totaloh']) %in% c( "integer", "numeric")){
      stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
    }
  }


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

  return(transform_table)
}
