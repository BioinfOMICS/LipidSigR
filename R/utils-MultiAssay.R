#' @title as_MultiAssayExperiment
#' @description Combines two pre-constructed SummarizedExperiment objects
#' (e.g. lipidomics and RNA-seq) into a single MultiAssayExperiment object
#' for downstream multi-omics integration analysis. The colData of both
#' SummarizedExperiment objects must share the same column structure
#' (produced by constructing \code{lipidSE} with \code{se_type="de_two"}
#' or \code{se_type="de_multiple"}, and \code{rnaSE} with
#' \code{se_type="rna"} using a matching \code{group_info} layout, both
#' via \code{\link{as_summarized_experiment}}).
#'
#' @param lipidSE A SummarizedExperiment object constructed by
#' \code{as_summarized_experiment()} with \code{se_type="de_two"} or
#' \code{se_type="de_multiple"} (the only two layouts \code{rnaSE}'s
#' \code{group_info} can match; see \code{rnaSE} below).
#' @param rnaSE A SummarizedExperiment object constructed by
#' \code{as_summarized_experiment(..., se_type="rna")}, whose
#' \code{group_info} must follow the de_two or de_multiple column layout
#' (\code{"ml"}/\code{"corr"} layouts are not supported for \code{"rna"})
#' and must match the one used for \code{lipidSE}.
#'
#' @return A MultiAssayExperiment object.
#' @export
#' @examples
#' library(dplyr)
#' library(magrittr)
#' data("abundance_lipidOmics")
#' data("abundance_rna")
#' data("group_info_multiOmics")
#'
#' ## lipidomics SE - de_two-style group_info (2 groups)
#' parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_lipidOmics$feature)
#' recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
#' lipid_abundance <- abundance_lipidOmics %>% dplyr::filter(feature %in% recognized_lipid)
#' lipid_char_table <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipid)
#' lipidSE <- as_summarized_experiment(
#'     lipid_abundance, lipid_char_table, group_info=group_info_multiOmics,
#'     se_type='de_two', paired_sample=FALSE)
#'
#' ## RNA SE built with se_type='rna', reusing the same group_info so its
#' ## colData structure matches lipidSE's
#' rnaSE <- as_summarized_experiment(
#'     abundance_rna, goslin_annotation=NULL, group_info=group_info_multiOmics,
#'     se_type='rna')
#'
#' ## combine into a single MultiAssayExperiment
#' mae <- as_MultiAssayExperiment(lipidSE, rnaSE)
as_MultiAssayExperiment <- function(lipidSE, rnaSE) {

    ## --- check both inputs are SummarizedExperiment objects ---
    .check_inputSE(lipidSE, metadata_list=NULL)
    .check_inputSE(rnaSE, metadata_list=NULL)

    ## --- check sample name overlap between the two experiments ---
    common_samples <- intersect(colnames(lipidSE), colnames(rnaSE))
    if (length(common_samples) == 0) {
        stop("No common sample names found between lipidSE and rnaSE. Please check the colnames of both objects.")
    }
    if (!setequal(colnames(lipidSE), colnames(rnaSE))) {
        warning(
            "Samples are not identical between the two experiments. ",
            "Number of common samples: ", length(common_samples),
            ". The MultiAssayExperiment will still be constructed, but ",
            "non-overlapping samples will appear as NA in the corresponding assay.")
    }

    ## --- validate and construct subject-level colData ---
    mae_colData <- .check_mae_colData(lipidSE, rnaSE, common_samples)

    ## --- construct ExperimentList ---
    expList <- MultiAssayExperiment::ExperimentList(
        list(lipidomics=lipidSE, rna=rnaSE))

    ## --- construct MultiAssayExperiment ---
    ## sampleMap is auto-generated when colnames match colData rownames
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments=expList,
        colData=mae_colData)

    ## --- summary message ---
    message(.mae_info(mae))

    return(mae)
}

#' @noRd
#' @description Validates that colData structures from lipidSE and rnaSE
#' are consistent (same columns, same values for common samples), and
#' returns a subject-level data frame covering the union of samples from
#' both experiments - not just the common ones - so that
#' MultiAssayExperiment's automatic sampleMap generation can find every
#' assay column name among the colData rownames (samples unique to one
#' experiment then correctly surface as NA for the other experiment).
.check_mae_colData <- function(lipidSE, rnaSE, common_samples) {

    ## --- colData must carry real group_info-derived columns (i.e. not
    ## built with se_type="profiling"). "profiling" stores colnames as a
    ## plain character vector, which SummarizedExperiment() still coerces
    ## into a 1-column colData (later force-renamed to "sample_name" by
    ## .extract_df() below) - so this is a "<= 1" check, not "== 0". Every
    ## valid de_two/de_multiple layout has at least 3 columns
    ## (sample_name, label_name, group[, pair]). ---
    lipid_colData <- SummarizedExperiment::colData(lipidSE)
    rna_colData <- SummarizedExperiment::colData(rnaSE)

    if (ncol(lipid_colData) <= 1) {
        stop(
            "lipidSE has no group_info-derived colData columns - it looks ",
            "like it was built with se_type='profiling'. ",
            "as_MultiAssayExperiment() requires lipidSE to be built with ",
            "se_type='de_two' or se_type='de_multiple'.")
    }
    if (ncol(rna_colData) <= 1) {
        stop(
            "rnaSE has no group_info-derived colData columns. ",
            "as_MultiAssayExperiment() requires rnaSE to be built with ",
            "se_type='rna' and a de_two or de_multiple-style group_info.")
    }

    ## --- colData must have a real 'sample_name' column. This has to be
    ## checked on the RAW colData, before .extract_df() below - .extract_df()
    ## unconditionally force-renames the first column to "sample_name" (it
    ## assumes the SE was built by as_summarized_experiment(), where that's
    ## always true), so checking colnames(lipid_group)/colnames(rna_group)
    ## afterwards would never actually catch a missing 'sample_name' column. ---
    if (isFALSE("sample_name" %in% colnames(as.data.frame(lipid_colData)))) {
        stop("colData of lipidSE must contain a 'sample_name' column.")
    }
    if (isFALSE("sample_name" %in% colnames(as.data.frame(rna_colData)))) {
        stop("colData of rnaSE must contain a 'sample_name' column.")
    }

    lipid_group <- .extract_df(lipidSE, type="group")
    rna_group <- .extract_df(rnaSE, type="group")

    ## --- colData column structure must match between the two SEs ---
    ## a mismatch here typically means lipidSE and rnaSE were constructed
    ## with incompatible se_type values (e.g. "de_two" vs "ml")
    if (isFALSE(identical(sort(colnames(lipid_group)), sort(colnames(rna_group))))) {
        stop(
            "colData column structure differs between lipidSE and rnaSE ",
            "(lipidSE: ", paste(colnames(lipid_group), collapse=", "),
            "; rnaSE: ", paste(colnames(rna_group), collapse=", "),
            "). This usually means the two objects were constructed with ",
            "incompatible se_type values. Please ensure both ",
            "SummarizedExperiment objects are constructed with a ",
            "compatible se_type so their colData structures match.")
    }

    ## --- align column order so lipid/rna group tables are comparable ---
    rna_group <- rna_group %>% dplyr::select(dplyr::all_of(colnames(lipid_group)))

    lipid_common <- lipid_group %>%
        dplyr::filter(sample_name %in% common_samples) %>%
        dplyr::arrange(sample_name)
    rna_common <- rna_group %>%
        dplyr::filter(sample_name %in% common_samples) %>%
        dplyr::arrange(sample_name)

    ## --- values must be consistent for common samples ---
    if (isFALSE(isTRUE(all.equal(lipid_common, rna_common, check.attributes=FALSE)))) {
        stop(
            "colData values differ between lipidSE and rnaSE for common ",
            "samples. Please ensure both SummarizedExperiment objects were ",
            "constructed from the same group_info table.")
    }

    ## --- build colData over the UNION of samples, not just the common
    ## ones. lipid_group already covers every lipidSE sample (including all
    ## common ones, whose values were just verified to match rnaSE); we
    ## only need to append the rows for samples that exist in rnaSE but not
    ## in lipidSE. This ensures every column name of both experiments is
    ## present in colData's rownames, which MultiAssayExperiment's implicit
    ## sampleMap generation requires - without this, samples unique to one
    ## experiment would be silently dropped instead of surfacing as NA. ---
    rna_only <- rna_group %>%
        dplyr::filter(!(sample_name %in% lipid_group$sample_name))

    mae_colData <- dplyr::bind_rows(lipid_group, rna_only) %>%
        dplyr::arrange(sample_name) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var="sample_name")

    return(mae_colData)
}

.mae_info <- function(mae) {
    note <- "MultiAssayExperiment info \n"
    note <- append(note, paste0(
        "Experiments: ", paste(names(mae), collapse=", "), "\n"))
    for (exp_name in names(mae)) {
        note <- append(note, paste0(
            "  ", exp_name, ": ",
            nrow(mae[[exp_name]]), " features x ",
            ncol(mae[[exp_name]]), " samples\n"))
    }
    note <- append(note, paste0(
        "Number of subjects with complete data across all experiments: ",
        sum(MultiAssayExperiment::complete.cases(mae))))
    return(note)
}
