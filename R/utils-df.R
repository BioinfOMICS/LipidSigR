#' @title as_summarized_experiment
#' @description Constructs a \code{SummarizedExperiment} object from input
#' data frames for downstream lipidomics analysis.
#'
#' @param abundance A data frame of lipid abundance data with the following
#' structure:
#' \itemize{
#'   \item \strong{Column 1} (\code{feature}): a character column of unique
#'   lipid species names used for feature-level analysis and lipid
#'   recognition/mapping. Must be named \code{"feature"}. No duplicate or
#'   missing values are allowed.
#'   \item \strong{Columns 2 to N}: numeric columns of sample abundance values,
#'   one column per sample. Column names must match the \code{sample_name}
#'   values in \code{group_info}. \code{NA}s are allowed.
#' }
#'
#' @param goslin_annotation A data frame of lipid structural annotations
#' returned by \code{\link[rgoslin]{parseLipidNames}} from the
#' rgoslin package. Before passing this data frame to the
#' function, remove any rows where the \code{Grammar} column is
#' \code{'NOT_PARSEABLE'}, as these entries represent lipid names that could
#' not be mapped to a recognised lipid grammar and are not supported for
#' downstream analysis. The corresponding lipids should also be removed from
#' \code{abundance}. Required for all lipidomics \code{se_type} options
#' (\code{"profiling"}, \code{"de_two"}, \code{"de_multiple"}, \code{"ml"},
#' \code{"corr"}). Must be \code{NULL} when \code{se_type = "rna"}, since
#' non-lipidomics abundance data has no lipid structural annotation.
#'
#' @param group_info A data frame containing sample metadata. The required
#' structure depends on the \code{se_type} selected:
#' \describe{
#'   \item{\code{"profiling"}}{No \code{group_info} is required; this argument
#'   is optional and can be omitted.}
#'
#'   \item{\code{"de_two"}}{Required. Compares lipid abundance between exactly
#'   two groups, with optional paired design support. Must contain the
#'   following columns in order: \code{sample_name}, \code{label_name},
#'   \code{group}, \code{pair}.
#'   \itemize{
#'     \item \code{sample_name}: unique sample identifiers that match the
#'     column names in \code{abundance} exactly. No \code{NA}s allowed.
#'     \item \code{label_name}: human-readable label for each sample
#'     (e.g., \code{ctrl1}, \code{hfref9}). No \code{NA}s allowed.
#'     \item \code{group}: biological group used for differential abundance
#'     testing. Must contain exactly 2 unique values. No \code{NA}s allowed.
#'     \item \code{pair}: paired design indicator. If samples are paired,
#'     provide sequential integer identifiers (e.g., 1 to N) with no gaps or
#'     missing values. If samples are unpaired, set all values to \code{NA}.
#'   }}
#'
#'   \item{\code{"de_multiple"}}{Required. Same structure as \code{"de_two"}
#'   but for comparisons across more than two groups. The \code{group} column
#'   must contain more than 2 unique values. The \code{pair} column is not
#'   required and should be omitted. Must contain the following columns in
#'   order: \code{sample_name}, \code{label_name}, \code{group}.}
#'
#'   \item{\code{"ml"}}{Required. Provides the outcome label for model
#'   training. Must contain the following columns in order:
#'   \code{sample_name}, \code{group}.
#'   \itemize{
#'     \item \code{sample_name}: unique sample identifiers that match the
#'     column names in \code{abundance} exactly. No \code{NA}s allowed.
#'     \item \code{group}: the class label the model learns to predict.
#'     For binary classification, must be encoded as \code{0}/\code{1}.
#'     No \code{NA}s allowed. Note: this is a prediction outcome label,
#'     not an experimental group for statistical testing.
#'   }}
#'
#'   \item{\code{"corr"}}{Required. Associates lipid abundance with continuous
#'   or categorical phenotype variables for correlation and/or regression
#'   analysis. Must contain the following columns:
#'   \itemize{
#'     \item \code{sample_name}: unique sample identifiers that match the
#'     column names in \code{abundance} exactly. No \code{NA}s allowed.
#'     \item Two or more phenotype or clinical variable columns (e.g.,
#'     \code{Age}, \code{BMI}, \code{Smoking}). Unlike other \code{se_type}
#'     options, there is no single \code{group} column; multiple phenotype
#'     columns are provided instead. Categorical variables must be pre-encoded
#'     as numeric values (e.g., binary categories as \code{0}/\code{1})
#'     before input.
#'   }}
#'
#'   \item{\code{"rna"}}{Required. \code{se_type = "rna"} builds a
#'   \code{SummarizedExperiment} for non-lipidomics abundance data (e.g. RNA
#'   expression) that is meant to feed downstream Multi-Omics analysis (see
#'   \code{\link{as_MultiAssayExperiment}}) rather than one of the 4
#'   lipidomics analysis designs above. \code{group_info} must be provided;
#'   its format is auto-detected from the column layout and must match one
#'   of the following two lipidomics layouts (\code{"ml"}/\code{"corr"}
#'   layouts are not supported for \code{"rna"}):
#'   \itemize{
#'     \item 4 columns, \code{sample_name}, \code{label_name}, \code{group},
#'     \code{pair} (same layout as \code{"de_two"}): \code{group} must
#'     contain exactly 2 unique values; paired design is inferred directly
#'     from the \code{pair} column (all \code{NA} = unpaired, sequential
#'     integers = paired) rather than a separate \code{paired_sample}
#'     argument.
#'     \item 3 columns, \code{sample_name}, \code{label_name}, \code{group}
#'     (same layout as \code{"de_multiple"}): \code{group} must contain more
#'     than 2 unique values, each with more than 2 samples.
#'   }
#'   In both cases \code{sample_name} must match the column names in
#'   \code{abundance} exactly with no duplicates, and no \code{NA}s are
#'   allowed in \code{sample_name}, \code{label_name}, or \code{group}.}
#' }
#'
#' @param se_type A character string specifying the analysis type for which
#' the \code{SummarizedExperiment} object is constructed. Default is \code{'profiling'}.
#' Must be one of:
#' \itemize{
#'   \item \code{"profiling"}: lipid profiling and visualisation.
#'   \item \code{"de_two"}: differential abundance analysis between two groups.
#'   \item \code{"de_multiple"}: differential abundance analysis across
#'   multiple groups.
#'   \item \code{"ml"}: machine learning classification.
#'   \item \code{"corr"}: correlation and regression analysis with phenotype
#'   variables.
#'   \item \code{"rna"}: non-lipidomics (e.g. RNA) abundance data intended for
#'   downstream Multi-Omics analysis. When \code{se_type = "rna"},
#'   \code{goslin_annotation} is not used and must be \code{NULL} (no lipid
#'   characteristics table is built).
#' }
#'
#' @param paired_sample Logical or \code{NULL}. For two-group analysis
#' (\code{se_type = "de_two"}), set to \code{TRUE} if the samples are paired
#' or \code{FALSE} if unpaired. For all other \code{se_type} options
#' (including \code{"de_multiple"}, \code{"ml"}, \code{"corr"}), set to
#' \code{NULL}. For \code{se_type = "rna"} the argument is not used and any
#' value is ignored.
#'
#' @importFrom rgoslin parseLipidNames
#' @return A \code{SummarizedExperiment} object configured for the specified
#' analysis type.
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
        abundance, goslin_annotation=NULL, group_info=NULL, se_type='profiling',
        paired_sample=FALSE){
    ## check input data --------------
    match.arg(se_type, c("profiling", "de_two", "de_multiple", "ml", "corr", 'rna'))
    .check_valid_inputData(
        abundance, goslin_annotation, group_info, se_type, paired_sample)

    abundance <- as.data.frame(abundance)

    ## ---- non-lipidomics path (e.g. se_type='rna') ----
    ## No goslin_annotation / lipid characteristics table exists for these
    ## omics, so skip all lipid-specific annotation/matching steps. The
    ## resulting SE feeds into Multi-Omics analysis rather than the 5
    ## lipidomics analysis designs above, so group_info is optional and only
    ## checked for sample alignment.
    if (se_type == "rna") {
        return(.as_summarized_experiment_nonLipid(abundance, group_info, se_type))
    }

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
            stop('Please verify if the `abundance` and the `goslin_annotation` ',
                 'tables originate from the same dataset, as all lipid names ',
                 'appear to be different between the two.')
        }
        if(isFALSE(isTRUE(all.equal(goslin_annotation$Original.Name, con.feature)))) {
            goslin_diff_lipid <- length(setdiff(
                goslin_annotation$Original.Name, con.feature))
            goslin_annotation %<>% dplyr::filter(Original.Name %in% con.feature)
            warning(goslin_diff_lipid, ' Original.Names have been removed from ',
                    'the `goslin_annotation` table as they differed from the ',
                    'features in the `abundance` table.')
        }
        if(isFALSE(isTRUE(all.equal(abundance_processed$feature, con.feature)))){
            abundance_diff_lipid <- length(setdiff(
                abundance_processed$feature, con.feature))
            abundance_processed %<>% dplyr::filter(feature %in% con.feature)
            warning(abundance_diff_lipid, ' features have been removed from ',
                    'the `abundance` table as they differed from the ',
                    'Original.Names in the `goslin_annotation` table.')
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

## Constructs a SummarizedExperiment for non-lipidomics abundance data (e.g.
## se_type='rna'), where no goslin_annotation / lipid characteristics table
## is available. rowData only carries the feature names; group_info is
## required and must follow the de_two or de_multiple column layout (checked
## by .check_rna_group_info()) - it is used to align/order samples in
## colData. No downstream analysis design is enforced beyond that, since the
## SE is meant to feed downstream Multi-Omics analysis rather than one of
## the 5 lipidomics analysis designs.
.as_summarized_experiment_nonLipid <- function(abundance, group_info, se_type) {
    ## check abundance data (feature name / numeric / duplicate / constant checks)
    abundance_processed <- .check_abundance(abundance, se_type)

    ## check group_info (required) and set colData
    group_info <- as.data.frame(group_info)
    detected_format <- .check_rna_group_info(abundance_processed, group_info)
    sample_order <- colnames(abundance_processed[-1])
    seCol <- group_info[match(sample_order, group_info$sample_name), , drop=FALSE]

    ## final check for sample number (de_two and de_multiple both require
    ## at least 2 samples, same as the matching lipidomics se_type)
    sample_limit <- 2
    if (isFALSE(.check_colNum(abundance_processed[-1], sample_limit)) ) {
        stop("Not enough samples. At least ", sample_limit, " samples are required.")
    }

    ## construct SE objects
    abundance_mat <- abundance_processed %>% dplyr::arrange(feature) %>%
        tibble::column_to_rownames(var="feature")
    feature_char <- data.frame(feature=rownames(abundance_mat))

    SE_data <- SummarizedExperiment::SummarizedExperiment(
        assays=list(abundance=as.matrix(abundance_mat)),
        rowData=S4Vectors::DataFrame(feature_char, row.names=feature_char$feature),
        colData=seCol)
    ## summary
    message(.inputData_info(SE_data, group_info, se_type, paired_sample=NULL))
    return(SE_data)
}

.inputData_info <- function(se, group_info, se_type, paired_sample) {
    note <- "Input data info \n"
    note <- append(note, paste0("se_type: ", se_type, "\n"))
    feature_label <- if (se_type=="rna") "features" else "lipids (features)"
    note <- append(note, paste0(
        "Number of ", feature_label, " available for analysis: ", dim(se)[1], "\n"))
    note <- append(note, paste0("Number of samples: ", dim(se)[2], "\n"))
    if (se_type=="rna") {
        ## 'rna' group_info follows either the de_two or de_multiple column
        ## layout (validated by .check_rna_group_info()), so "group" and,
        ## when present, "pair" can be reported the same way.
        if (!is.null(group_info)) {
            se_group_info <- as.data.frame(SummarizedExperiment::colData(se))
            note <- append(note, paste0("Number of group: ", length(unique(se_group_info$group)), "\n"))
            if ("pair" %in% colnames(se_group_info)) {
                if (isFALSE(isTRUE(all(is.na(se_group_info$pair))))) {
                    paired <- "Paired samples."
                } else {
                    paired <- "Not paired samples."
                }
                note <- append(note, paired)
            }
        } else {
            note <- append(note, "Number of group: Group information table is not provided.")
        }
    } else if (se_type!="corr") {
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
#' LipidSigR, such as Profiling, Enrichment and subChar.
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
        dplyr::group_by(char) %>%
        dplyr::filter(dplyr::n_distinct(char_feature) > 1) %>%
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
