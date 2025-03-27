#' @title as_summarized_experiment
#' @description This function construct a SummarizedExperiment object for analysis.
#' by input data frames.
#' @param abundance A data frame of predictors includes features (such
#' as molecules or lipid classes) and their abundance in each sample. NAs are
#' allowed in all columns except the first. The first column must be named "feature"
#' and contain names of unique lipid species.
#' @param goslin_annotation A data frame of lipid characteristics from the
#' \code{\link[rgoslin]{parseLipidNames}} function of \pkg{\link{rgoslin}} package.
#' Please exclude lipids with the "Grammar" column marked as 'NOT_PARSEABLE'.
#' @param group_info A data frame containing group information. The required
#' columns are listed as below.
#' \enumerate{
#' \item sample_name: name of each sample. NAs are not allowed; the sample names
#' must be unique and match those in the abundance data.
#' \item label_name: label name of each sample. NAs are not allowed.
#' \item group: group name of each sample. NAs are not allowed.
#' \item pair: (only required for two-group data) If paired samples are contained,
#' they should be sequentially numbered from 1 to N without missing, blank, or
#' skipped numbers. If not, the values must all be marked as NA.
#' }
#' @param se_type Character. The SummarizedExperiment object is generated for
#' running functions related to which analysis type. Must by one of "profiling",
#' "de_two", "de_multiple", "ml", or "corr".
#' @param paired_sample Logical/NULL. For two-group data, enter TRUE/FALSE
#' to indicate if paired samples are included. For multi-group data, enter NULL.
#' @importFrom rgoslin parseLipidNames
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' library(dplyr)
#' data("abundance_twoGroup")
#' data("group_info_twoGroup")
#' parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_twoGroup$feature)
#' recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
#' abundance <- abundance_twoGroup %>% dplyr::filter(feature %in% recognized_lipid)
#' goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipid)
#' se <- as_summarized_experiment(abundance, goslin_annotation,
#'      group_info=group_info_twoGroup, se_type='de_two', paired_sample=FALSE)
as_summarized_experiment <- function(
        abundance, goslin_annotation, group_info=NULL,
        se_type=c('profiling', 'de_two', 'de_multiple', 'ml', 'corr'),
        paired_sample=FALSE){
    ## check input data --------------
    if (is.null(se_type) | isFALSE(se_type %in% c("profiling", "de_two", "de_multiple", "ml", "corr")) ) {
        stop("se_type must be one of 'profiling', 'de_two', 'de_multiple', 'ml', or 'corr'.")
    }
    .check_valid_inputData(
        abundance, goslin_annotation, group_info, se_type, paired_sample)

    abundance <- as.data.frame(abundance)
    goslin_annotation <- as.data.frame(goslin_annotation)

    ## check goslin_annotation not duplicate & table column
    goslin_annotation <- .check_goslin_annotation(goslin_annotation)

    ## check group_info
    if (se_type %in% c("de_two", "de_multiple", "ml", "corr") ) {
        group_info <- as.data.frame(group_info)
        .check_group_info(abundance, group_info, se_type, paired_sample)
    }
    ## check abundance data
    abundance_processed <- .check_abundance(abundance, se_type)
    ## Goslin lipid char table match Abundance data
    if (isFALSE(isTRUE(all.equal(
        goslin_annotation$Original.Name, abundance_processed$feature)))) {
        con.feature <- intersect(
            goslin_annotation$Original.Name, abundance_processed$feature)
        if(length(con.feature) == 0){
            stop(paste0(
                'Please verify if the `abundance` and the `goslin_annotation` ',
                'tables originate from the same dataset, as all lipid names ',
                'appear to be different between the two.'))
        }
        if(isFALSE(isTRUE(all.equal(goslin_annotation$Original.Name, con.feature)))) {
            goslin_diff_lipid <- length(setdiff(
                goslin_annotation$Original.Name, con.feature))
            goslin_annotation %<>% dplyr::filter(Original.Name %in% con.feature)
            warning(paste0(
                goslin_diff_lipid, ' Original.Names have been removed from ',
                'the `goslin_annotation` table as they differed from the ',
                'features in the `abundance` table.'))
        }
        if(isFALSE(isTRUE(all.equal(abundance_processed$feature, con.feature)))){
            abundance_diff_lipid <- length(setdiff(
                abundance_processed$feature, con.feature))
            abundance_processed %<>% dplyr::filter(feature %in% con.feature)
            warning(paste0(
                abundance_diff_lipid, ' features have been removed from ',
                'the `abundance` table as they differed from the ',
                'Original.Names in the `goslin_annotation` table.'))
        }
    }
    ## lipid annotation
    lipid_char <- lipid_annotation(goslin_annotation)
    ## Abundance data match group_info and set colData
    if (se_type!="profiling") {
        sample_order <- colnames(abundance[-1])
        seCol <- group_info[match(sample_order, group_info$sample_name), ]
    } else { ## Profiling
        seCol <- colnames(abundance)[-1]
    }
    ## final check for sample number
    sample_limit <- switch(
        se_type,
        "profiling"=2, "de_two"=2, "de_multiple"=2, "ml"=60, "corr"=10)
    if (isFALSE(.check_colNum(abundance[-1], sample_limit)) ) {
        stop("Not enough samples. Each se_type must at least have the following samples: profiling: 2, de_two: 2, de_multiple: 2, ml: 60, corr: 10.")
    }
    ## construct SE objects
    abundance_mat <- abundance_processed %>% dplyr::arrange(feature) %>%
        tibble::column_to_rownames(var="feature")
    lipid_char %<>% dplyr::arrange(feature)

    SE_data <- SummarizedExperiment::SummarizedExperiment(
        assays=list(abundance=as.matrix(abundance_mat) ),
        rowData=S4Vectors::DataFrame(lipid_char, row.names=lipid_char$feature),
        colData=seCol)
    ## summary
    message(.inputData_info(SE_data, group_info, se_type, paired_sample))
    return(SE_data)
}

.inputData_info <- function(se, group_info, se_type, paired_sample) {
    note <- "Input data info \n"
    note <- append(note, paste0("se_type: ", se_type, "\n"))
    note <- append(note, paste0("Number of lipids (features) available for analysis: ", dim(se)[1], "\n"))
    note <- append(note, paste0("Number of samples: ", dim(se)[2], "\n"))
    if (se_type!="corr") {
        if (!is.null(group_info)) {
            se_group_info <- as.data.frame(SummarizedExperiment::colData(se))
            note <- append(note, paste0("Number of group: ", length(unique(se_group_info$group)), "\n"))
            if (isTRUE(paired_sample) ){
                paired <- "Paired samples."
            } else {
                paired <- "Not paired samples."
            }
            note <- append(note, paired)
        } else {
            note <- append(note, "Number of group: Group information table is not provided.")
        }
    }
    return(note)
}

# description: get data.frames from SE
# input: 1 SE, type in (abundance, lipid, group, condition, adjusted)
# output: 1 data.frame
.extract_df <- function(se, type="abundance") {
    stopifnot(inherits(se, "SummarizedExperiment"))
    if (type=="abundance") {
        stopifnot("assays" %in% slotNames(se) )
        data <- as.data.frame(SummarizedExperiment::assay(se)) %>%
            tibble::rownames_to_column("feature")
    } else if (type=="lipid") {
        data <- as.data.frame(SummarizedExperiment::rowData(se))
    } else if (type=="group") {
        data <- as.data.frame(SummarizedExperiment::colData(se))
        colnames(data)[1] <- "sample_name"
    } else {
        stop("type should be one of abundance, lipid, or group")
    }
    return(data)
}

#' @title extract_summarized_experiment
#' @description This function returns a list of data frames includes abundance data, lipid
#' characteristic table, group information table, and other analysis results from
#' a SummarizedExperiment object.
#' @param se A SummarizedExperiment object.
#' @return Return a list of data frames.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
#'     significant='padj', p_cutoff=0.05, FC_cutoff=2, transform='log10')
#' table_list <- extract_summarized_experiment(deSp_se)
extract_summarized_experiment <- function(se){
    ## check
    if (isFALSE(inherits(se, "SummarizedExperiment") )) {
        stop("The input must be a SummarizedExperiment object construct
             by as_summarized_experiment function or ouput from upstream analysis function.")
    }
    if (all(dim(se) == c(0, 0) )) {
        stop("The input must be a SummarizedExperiment object construct by
             as_summarized_experiment function or ouput from upstream analysis function.")
    }
    if (any(
        length(SummarizedExperiment::assays(se))==0,
        length(SummarizedExperiment::rowData(se))==0,
        length(SummarizedExperiment::colData(se))==0) ) {
        stop("The input must be a SummarizedExperiment object construct by
             as_summarized_experiment function or ouput from upstream analysis function.")
    }
    ## correct slot
    slot_name <- c("assays", "colData", "elementMetadata")
    if (isFALSE(all(slot_name %in% slotNames(se)) )) {
        stop("Incorrect SummarizedExperiment structure, please construct by
             as_summarized_experiment function or ouput from upstream analysis function.")
    }
    abundance <- .extract_df(se, type = "abundance")
    lipid <- .extract_df(se, type = "lipid")
    group <- .extract_df(se, type = "group")
    metadata_se <- S4Vectors::metadata(se)

    if (ncol(group) > 1) {
        result_list <- c(
            list(abundance=abundance, lipid_char_table=lipid, group_info=group),
            metadata_se)
    } else {
        result_list <- c(
            list(abundance=abundance, lipid_char_table=lipid), metadata_se)
    }
    return(result_list=result_list)
}

#' @title list_lipid_char
#' @description This function lists the lipid characteristics suitable for
#' different types of analyses.
#' @param processed_se A SummarizedExperiment object after conducting \code{\link{data_process}}.
#' @return Return a list of named vector.
#' \enumerate{
#' \item deChar_list: List of lipid characteristics for conducting differential
#' expressed lipid characteristics analysis.
#' \item chain_db_list: List of lipid characteristics for performing chain
#' length-double bond correlation analysis.
#' \item common_list: List of lipid characteristics for other analyses in
#' LipidSigR, such as Profiling and Enrichment.
#' \item ml_char_list: List of lipid characteristics for conducting machine learning analysis.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' char_list <- list_lipid_char(processed_se)
list_lipid_char <- function(processed_se){
    #data("lipidChar")
    .check_inputSE(processed_se, metadata_list=NULL)
    ## read files
    abundance <- .extract_df(processed_se, type='abundance') %>%
        tidyr::gather(sample_name, values, -feature)
    lipid_char <- .extract_df(processed_se, type='lipid') %>%
        dplyr::select(feature, dplyr::all_of(lipidChar$characteristic)) %>%
        tidyr::gather(char, char_feature, -feature) %>%
        dplyr::filter(!is.na(char_feature))
    group_info <- .extract_df(processed_se, type='group')

    ## data-specific lipid characteristics
    data_info <- lipid_char %>%
        dplyr::left_join(abundance, by='feature', relationship='many-to-many') %>%
        dplyr::filter(!is.na(char_feature)) %>%
        dplyr::left_join(group_info, by='sample_name')

    ## for general situation
    common_char <- data_info %>%
        dplyr::group_by(char, char_feature, sample_name) %>%
        dplyr::summarise(sum(values, na.rm=TRUE), .groups='drop') %>%
        dplyr::distinct(char) %>%
        dplyr::left_join(lipidChar, by=c('char'='characteristic')) %>%
        dplyr::arrange(match(
            aspect, c('Lipid classification', 'Fatty acid properties',
                      'Physical or chemical properties', 'Cellular component',
                      'Function')), char)
    common_char_list <- common_char$char
    names(common_char_list) <- common_char$aspect

    if (ncol(group_info) > 1) {
        if (all(
            suppressWarnings(
                unique(group_info$group)==c(0,1)
            )
            )) {
            ml_char_list <- common_char_list[common_char_list %in% c("class", "Total.OH", "Total.DB", "FA.OH", "FA.DB", "FA.C")]
            de_char_list <- "Only provided when differential analysis data is input."
            heatmap_char_list <- "Only provided when differential analysis data is input."
        } else {
            ml_char_list <- "Only provided when machine learning analysis data is input."
            ## data-specific ratio characteristics
            sp2ratio <- convert_sp2ratio(processed_se, transform='none')
            data_ratio <- S4Vectors::metadata(sp2ratio)[["char_list"]]
            names(data_ratio) <- rep('Specific ratios', length(data_ratio))
            ## for DE char
            de_char_list <- c(common_char_list, data_ratio)
            ## for heatmap_chain_db
            heatmap_char <- common_char %>%
                dplyr::filter(aspect != 'Fatty acid properties')
            heatmap_char_list <- heatmap_char$char
            names(heatmap_char_list) <- heatmap_char$aspect
        }
    } else {
        de_char_list <- "Only provided when differential analysis data is input."
        heatmap_char_list <- "Only provided when differential analysis data is input."
        ml_char_list <- "Only provided when machine learning analysis data is input."
    }
    return(list(
        deChar_list=de_char_list, chain_db_list=heatmap_char_list,
        common_list=common_char_list, ml_char_list=ml_char_list))
}
