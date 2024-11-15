#' @title convert_sp2char
#' @description This function sums up lipid species abundance according to lipid characteristics,
#' which converts the lipid species abundance table to the lipid characteristic abundance table.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", and "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' char_data <- convert_sp2char(processed_se, transform='log10')

convert_sp2char <- function(
        processed_se, transform=c('none', 'log10', 'square', 'cube')){

    # check input SE
    .check_inputSE(processed_se, metadata_list=NULL)
    abundance <- .extract_df(processed_se, type = "abundance")
    abundance[-1] <- lapply(abundance[-1], as.numeric)
    lipid_char_table <- .extract_df(processed_se, type = "lipid")
    group_info <- .extract_df(processed_se, type = "group")

    char_var_list <- lipidChar$characteristic

    .check_imputation(abundance)
    if (isTRUE(.check_nor_negative(abundance[-1])) ) {
      stop("Detect negative values in the input abundance data. We recommended choosing a different normalization method during data processing.")
    }

    # check parameter
    if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
      stop("transform must be one of 'none', 'log10', 'square' or 'cube'.")
    }

    ## sum each char exp
    char_data_all <- data.frame()
    for (i in seq(char_var_list)) {
        char <- char_var_list[i]
        ## sum characteristics
        abundance_ga <- abundance %>%
            dplyr::left_join(
                lipid_char_table[c('feature', char)], by='feature') %>%
            dplyr::select(-1)
        if (nrow(abundance_ga) == 0){
            char_data <- NULL
        } else if (all(is.na(abundance_ga[, char]))) {
            char_data <- NULL
        } else {
            char_data <- .sum_by_char(
                abundance_ga, lipid_char_table, char, transform)
            char_data <- char_data %>% dplyr::rename(feature=char) %>%
                dplyr::mutate(characteristic = char, .before=feature)
        }
        ## binding long data frame
        if (i==1) {
            char_data_all <- char_data
        } else {
            char_data_all <- rbind(char_data_all, char_data)
        }
    }
    if(is.null(char_data_all)) {
        stop("No enough data.")
    }
    transform_data <- char_data_all %>%
        dplyr::mutate(rowName=paste(characteristic, feature, sep="|"))
    transform_select <- transform_data %>%
        tibble::column_to_rownames(var="rowName")
    transform_mat <- as.matrix(transform_select[, group_info$sample_name])
    row_data <- char_data_all[, c("characteristic", "feature")]

    char_se <- SummarizedExperiment::SummarizedExperiment(
        assays=list(transform_table=transform_mat),
        rowData=row_data,
        colData=SummarizedExperiment::colData(processed_se))

  return(char_se=char_se)
}

.sum_by_char <- function(abundance_ga, lipid_char_table, char, transform){
    ## process "|" characteristic
    if(isTRUE(any(stringr::str_detect(abundance_ga[, char], '\\|')) )){
        transform_table <- tidyr::separate_rows(
            abundance_ga, char, sep = "\\|")
    } else {
        transform_table <- abundance_ga
    }
    transform_table <- transform_table %>%
        stats::aggregate(
            stats::as.formula(stringr::str_c('. ~ ', char)), ., sum)
    ## data process: transformation
    transform_table_process <- .transform(transform_table, transform)
    return(char_data=transform_table_process)
}
