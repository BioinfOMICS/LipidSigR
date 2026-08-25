library(testthat)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(rgoslin)

## --- shared fixtures -------------------------------------------------
data("abundance_lipidOmics")
data("abundance_rna")
data("group_info_multiOmics")

parse_lipid <- rgoslin::parseLipidNames(lipidNames = abundance_lipidOmics$feature)
recognized_lipid <- parse_lipid$Original.Name[
    which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
lipid_abundance <- abundance_lipidOmics %>%
    dplyr::filter(feature %in% recognized_lipid)
lipid_char_table <- parse_lipid %>%
    dplyr::filter(Original.Name %in% recognized_lipid)

## lipidSE and rnaSE built from the same group_info (de_two format, 2
## groups: tsne0/tsne1) and the exact same 59 samples -> full overlap.
lipidSE <- as_summarized_experiment(
    lipid_abundance, lipid_char_table, group_info = group_info_multiOmics,
    se_type = 'de_two', paired_sample = FALSE)
rnaSE <- as_summarized_experiment(
    abundance_rna, goslin_annotation = NULL, group_info = group_info_multiOmics,
    se_type = 'rna')

expect_correct_mae <- function(mae) {
    expect_true(inherits(mae, "MultiAssayExperiment"))
    expect_setequal(names(mae), c("lipidomics", "rna"))
    expect_equal(nrow(mae[["lipidomics"]]), nrow(lipidSE))
    expect_equal(nrow(mae[["rna"]]), nrow(rnaSE))
}

## Hand-built minimal SummarizedExperiment (bypassing as_summarized_experiment())
## so we can hand it an arbitrary colData - used to reach the defensive
## colData checks in .check_mae_colData() that as_summarized_experiment()
## itself can never actually produce (it always yields either a proper
## sample_name-first colData, or the 1-column profiling one).
make_minimal_SE <- function(sample_names, colData_df) {
    n_sample <- length(sample_names)
    mat <- matrix(
        seq_len(2 * n_sample), nrow = 2,
        dimnames = list(c("f1", "f2"), sample_names))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = mat),
        rowData = S4Vectors::DataFrame(
            feature = c("f1", "f2"), row.names = c("f1", "f2")),
        colData = colData_df)
}

# Test for basic functionality and output structure -----------------------
test_that("as_MultiAssayExperiment builds a MultiAssayExperiment when samples fully overlap", {
    mae <- as_MultiAssayExperiment(lipidSE, rnaSE)
    expect_correct_mae(mae)
    ## full overlap -> colData should cover exactly the shared sample set
    expect_setequal(
        rownames(MultiAssayExperiment::colData(mae)), colnames(lipidSE))
    expect_setequal(
        rownames(MultiAssayExperiment::colData(mae)), colnames(rnaSE))
    ## every column of both experiments should be represented in sampleMap
    smap <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    expect_setequal(smap$colname[smap$assay == "lipidomics"], colnames(lipidSE))
    expect_setequal(smap$colname[smap$assay == "rna"], colnames(rnaSE))
    ## every subject has complete data across both experiments
    expect_equal(
        sum(MultiAssayExperiment::complete.cases(mae)),
        length(colnames(lipidSE)))
})

test_that("as_MultiAssayExperiment prints a summary message", {
    expect_message(
        as_MultiAssayExperiment(lipidSE, rnaSE), "MultiAssayExperiment info")
})

# Test for sample overlap handling -----------------------------------------
test_that("as_MultiAssayExperiment errors when there is no sample overlap", {
    synthetic_abundance <- data.frame(
        feature = c("geneA", "geneB", "geneC"),
        s1 = c(1, 2, 3), s2 = c(2, 3, 4), s3 = c(3, 4, 5), s4 = c(4, 5, 6))
    synthetic_group <- data.frame(
        sample_name = c("s1", "s2", "s3", "s4"),
        label_name = c("s1", "s2", "s3", "s4"),
        group = c("A", "A", "B", "B"),
        pair = NA)
    rnaSE_noOverlap <- as_summarized_experiment(
        synthetic_abundance, goslin_annotation = NULL,
        group_info = synthetic_group, se_type = 'rna')

    expect_snapshot(
        as_MultiAssayExperiment(lipidSE, rnaSE_noOverlap), error=TRUE)
})

test_that("as_MultiAssayExperiment warns on partial overlap and covers the union of samples in colData", {
    keep_samples <- colnames(rnaSE)[1:50]
    rnaSE_partial <- rnaSE[, keep_samples]
    lipid_only <- setdiff(colnames(lipidSE), keep_samples)
    expect_true(length(lipid_only) > 0)

    expect_warning(
        mae <- as_MultiAssayExperiment(lipidSE, rnaSE_partial),
        "Samples are not identical between the two experiments")

    union_samples <- union(colnames(lipidSE), keep_samples)
    ## colData must cover the union, not just the intersection, of samples
    expect_setequal(
        rownames(MultiAssayExperiment::colData(mae)), union_samples)

    ## sampleMap must still list every column from both experiments, so
    ## lipidSE-only samples are not silently dropped
    smap <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    expect_setequal(smap$colname[smap$assay == "lipidomics"], colnames(lipidSE))
    expect_setequal(smap$colname[smap$assay == "rna"], keep_samples)

    ## lipid-only samples have no rna data -> incomplete cases
    expect_true(sum(MultiAssayExperiment::complete.cases(mae)) < length(union_samples))
})

# Test for colData structure/consistency validation ------------------------
test_that("as_MultiAssayExperiment errors when lipidSE has no group_info-derived colData (se_type = 'profiling')", {
    lipidSE_profiling <- as_summarized_experiment(
        lipid_abundance, lipid_char_table, group_info = NULL,
        se_type = 'profiling', paired_sample = NULL)
    expect_equal(ncol(SummarizedExperiment::colData(lipidSE_profiling)), 1)
    expect_snapshot(
        as_MultiAssayExperiment(lipidSE_profiling, rnaSE), error=TRUE)
})

test_that("as_MultiAssayExperiment errors when rnaSE has no group_info-derived colData", {
    all_samples <- colnames(lipidSE)
    rnaSE_bad <- make_minimal_SE(
        all_samples,
        S4Vectors::DataFrame(sample_name = all_samples, row.names = all_samples))
    expect_equal(ncol(SummarizedExperiment::colData(rnaSE_bad)), 1)
    expect_snapshot(
        as_MultiAssayExperiment(lipidSE, rnaSE_bad), error=TRUE)
})

test_that("as_MultiAssayExperiment errors when lipidSE colData has no 'sample_name' column", {
    all_samples <- colnames(lipidSE)
    lipidSE_noSampleName <- make_minimal_SE(
        all_samples,
        S4Vectors::DataFrame(
            foo = all_samples, bar = rep("x", length(all_samples)),
            row.names = all_samples))
    expect_snapshot(
        as_MultiAssayExperiment(lipidSE_noSampleName, rnaSE), error=TRUE)
})

test_that("as_MultiAssayExperiment errors when rnaSE colData has no 'sample_name' column", {
    all_samples <- colnames(lipidSE)
    rnaSE_noSampleName <- make_minimal_SE(
        all_samples,
        S4Vectors::DataFrame(
            foo = all_samples, bar = rep("x", length(all_samples)),
            row.names = all_samples))
    expect_snapshot(
        as_MultiAssayExperiment(lipidSE, rnaSE_noSampleName), error=TRUE)
})

test_that("as_MultiAssayExperiment errors when lipidSE and rnaSE group_info structures differ", {
    group_info_deMultiple <- group_info_multiOmics %>%
        dplyr::select(-pair) %>%
        dplyr::mutate(group = rep(c("g1", "g2", "g3"), length.out = dplyr::n()))
    rnaSE_deMultiple <- as_summarized_experiment(
        abundance_rna, goslin_annotation = NULL,
        group_info = group_info_deMultiple, se_type = 'rna')

    expect_snapshot(
        as_MultiAssayExperiment(lipidSE, rnaSE_deMultiple), error=TRUE)
})

test_that("as_MultiAssayExperiment errors when common-sample group_info values differ between lipidSE and rnaSE", {
    group_info_flipped <- group_info_multiOmics
    group_info_flipped$group <- ifelse(
        group_info_flipped$group == "tsne0", "tsne1", "tsne0")
    rnaSE_flipped <- as_summarized_experiment(
        abundance_rna, goslin_annotation = NULL,
        group_info = group_info_flipped, se_type = 'rna')

    expect_snapshot(
        as_MultiAssayExperiment(lipidSE, rnaSE_flipped), error=TRUE)
})

# Test for wrong input types ------------------------------------------------
test_that("as_MultiAssayExperiment can handle wrong input", {
    expect_snapshot(as_MultiAssayExperiment(data.frame(), rnaSE), error=TRUE)
    expect_snapshot(as_MultiAssayExperiment(lipidSE, data.frame()), error=TRUE)
    expect_snapshot(as_MultiAssayExperiment(lipidSE, matrix(1:4, 2, 2)), error=TRUE)
    expect_snapshot(as_MultiAssayExperiment(NULL, rnaSE), error=TRUE)
    expect_snapshot(as_MultiAssayExperiment(lipidSE, NULL), error=TRUE)
})
