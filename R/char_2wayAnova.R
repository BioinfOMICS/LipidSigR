#' @title char_2wayAnova
#' @description This function calculates two-way ANOVA for all lipid
#' characteristics.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by
#' \code{\link{data_process}}.
#' @param ratio_transform Character. Method for transform the ratio-based
#' abundance. Allowed methods include "none" and "log2". Select 'none' to skip
#' data transformation. Default is \code{'log2'}.
#' @param char_transform Character. Method for transform the lipid
#' characteristics-based abundance. Allowed methods include "none", "log10",
#' "square", and "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return 1 table with two-way ANOVA results.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se_twoGroup <- data_process(
#'     se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' twoWayAnova_twoGroup <- char_2wayAnova(
#'     processed_se_twoGroup, ratio_transform='log2', char_transform='log10')
#' data("se_multiGroup")
#' processed_se_multiGroup <- data_process(
#'     se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' twoWayAnova_multiGroup <- char_2wayAnova(
#'     processed_se_multiGroup, ratio_transform='log2', char_transform='log10')
char_2wayAnova <- function(processed_se, ratio_transform=c('none', 'log2'),
                           char_transform=c('none', 'log10', 'square', 'cube')){

    ## Check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if(is.null(ratio_transform) |
       isFALSE(ratio_transform %in% c('none', 'log2'))){
        stop("The 'ratio_transform' parameter must be either 'none' or 'log2'.")
    }
    if(is.null(char_transform) |
       isFALSE(char_transform %in% c('none', 'log10', 'square', 'cube'))){
        stop("The 'char_transform' parameter must be one of 'none', 'log10',
             'square' or 'cube'.")
    }
    ## convert_sp2char
    char_se <- convert_sp2char(processed_se, transform=char_transform)
    ## convert_sp2ratio
    ratio_se <- convert_sp2ratio(processed_se, transform=ratio_transform)
    ## bind 2 SE
    SummarizedExperiment::assayNames(char_se) <-
        SummarizedExperiment::assayNames(ratio_se)
    all_char_se <- SummarizedExperiment::rbind(char_se, ratio_se)

    #### Extract data from SE ####
    ## Lipid abundance
    abundance <- .extract_df(all_char_se, type="abundance")
    ## Group information
    group_info <- .extract_df(all_char_se, type="group")
    ## Check nSample > 1
    nSample <- group_info %>% dplyr::group_by(group) %>%
        dplyr::summarise(n=dplyr::n(), .groups='drop')
    if(any(nSample$n < 2)){
        stop("This function requires each group to have at least two samples.")
    }
    ## Lipid characteristics table
    lipid_char <- .extract_df(all_char_se, type="lipid")

    ## Separate the feature
    abund_sep <- abundance %>%
        tidyr::separate(feature, into=c('characteristic', 'feature'), sep='\\|')
    ## Two-way ANOVA with interaction effect
    test_res <- .two_way_anova(abundance=abund_sep, group_info, char=NULL)
    test_res %<>%
        dplyr::left_join(lipidChar[, c('aspect', 'characteristic')],
                         by='characteristic') %>%
        dplyr::mutate(aspect=ifelse(is.na(aspect), 'Specific ratios', aspect)) %>%
        dplyr::select(aspect, dplyr::everything()) %>%
        dplyr::arrange(
            match(aspect, c('Lipid classification', 'Fatty acid properties',
                  'Physical or chemical properties', 'Cellular component',
                  'Function', 'Specific ratios'))) %>%
        as.data.frame()
    return(table_stat=test_res)
}
