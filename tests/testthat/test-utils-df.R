library(testthat)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(rgoslin)

data("abundance_twoGroup")
data("group_info_twoGroup")
two_parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_twoGroup$feature)
two_recognized_lipid <- two_parse_lipid$Original.Name[which(two_parse_lipid$Grammar != 'NOT_PARSEABLE')]
two_abundance <- abundance_twoGroup %>%
    dplyr::filter(feature %in% two_recognized_lipid)
two_char_table <- two_parse_lipid %>%
    dplyr::filter(Original.Name %in% two_recognized_lipid)

data("abundance_multiGroup")
data("group_info_multiGroup")
multi_parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_multiGroup$feature)
multi_recognized_lipid <- multi_parse_lipid$Original.Name[which(multi_parse_lipid$Grammar != 'NOT_PARSEABLE')]
multi_abundance <- abundance_multiGroup %>%
    dplyr::filter(feature %in% multi_recognized_lipid)
multi_char_table <- multi_parse_lipid %>%
    dplyr::filter(Original.Name %in% multi_recognized_lipid)

data("ml_abundance")
data("condition_table")
ml_parse_lipid <- rgoslin::parseLipidNames(lipidNames=ml_abundance$feature)
ml_recognized_lipid <- ml_parse_lipid$Original.Name[which(ml_parse_lipid$Grammar != 'NOT_PARSEABLE')]
ml_abundance <- ml_abundance %>%
    dplyr::filter(feature %in% ml_recognized_lipid)
ml_char_table <- ml_parse_lipid %>%
    dplyr::filter(Original.Name %in% ml_recognized_lipid)

data("corr_abundance")
data("corr_group_info")
corr_parse_lipid <- rgoslin::parseLipidNames(lipidNames=corr_abundance$feature)
corr_recognized_lipid <- corr_parse_lipid$Original.Name[which(corr_parse_lipid$Grammar != 'NOT_PARSEABLE')]
corr_abundance <- corr_abundance %>%
    dplyr::filter(feature %in% corr_recognized_lipid)
corr_char_table <- corr_parse_lipid %>%
    dplyr::filter(Original.Name %in% corr_recognized_lipid)


expect_correct_se <- function(se){
    expect_true(inherits(se, "SummarizedExperiment"))
    expect_true("assays" %in% slotNames(se),
                info = "Assays slot is missing")
    expect_true("colData" %in% slotNames(se),
                info = "colData slot is missing")
    expect_true("elementMetadata" %in% slotNames(se),
                info = "elementMetadata slot is missing")
}


# Test for basic functionality and output structure
test_that("as_summarized_experiment build SE correctly.", {
    profiling_se <- as_summarized_experiment(
        two_abundance, two_char_table, group_info=NULL, se_type='profiling', paired_sample=NULL)
    expect_correct_se(profiling_se)
    ## two
    two_se <- as_summarized_experiment(
        two_abundance, two_char_table, group_info_twoGroup, se_type='de_two',
        paired_sample=FALSE)
    expect_correct_se(two_se)
    group_info_twoGroup <- group_info_twoGroup[1:22, ] %>% dplyr::select(-pair) %>%
        dplyr::mutate(pair=rep(1:11, each = 2))
    two_se <- as_summarized_experiment(
        two_abundance[, 1:23], two_char_table, group_info_twoGroup, se_type='de_two',
        paired_sample=TRUE)
    expect_correct_se(two_se)
    ## multiple
    multi_se <- suppressWarnings(
        as_summarized_experiment(
            multi_abundance, multi_char_table, group_info_multiGroup, se_type='de_multiple',
            paired_sample=NULL))
    expect_correct_se(multi_se)
    multi_se <- suppressWarnings(
        as_summarized_experiment(
            multi_abundance, multi_char_table, group_info_multiGroup, se_type='de_multiple',
            paired_sample=NULL))
    expect_correct_se(multi_se)
    # ML
    ml_se <- suppressWarnings(
        as_summarized_experiment(
            ml_abundance, ml_char_table, condition_table, se_type='ml', paired_sample=NULL))
    expect_correct_se(ml_se)
    # Corr
    corr_se <- suppressWarnings(
        as_summarized_experiment(
            corr_abundance, corr_char_table, corr_group_info, se_type='corr', paired_sample=NULL))
    expect_correct_se(corr_se)
})

# Test for error input
test_that("as_summarized_experiment can handle NULL or 0*0 data frame error input.", {
    expect_error(as_summarized_experiment(
        abundance=NULL, two_char_table, group_info=NULL, se_type='profiling',
        paired_sample=FALSE), "cannot be null.")
    expect_error(
        as_summarized_experiment(
            two_abundance, goslin_annotation=NULL, group_info=NULL, se_type='profiling', paired_sample=FALSE),
        "Goslin lipid characteristics table cannot be null.")
    expect_error(
        as_summarized_experiment(
            abundance=data.frame(), two_char_table, group_info=NULL, se_type='profiling', paired_sample=FALSE),
        "must be a data frame and cannot be data frame with 0 columns and 0 rows.")
    expect_error(
        as_summarized_experiment(
            two_abundance, goslin_annotation=data.frame(), group_info=NULL, se_type='de_two', paired_sample=FALSE),
        "must be a data frame and cannot be data frame with 0 columns and 0 rows.")
    expect_error(as_summarized_experiment(
        two_abundance, two_char_table, group_info=data.frame(), se_type='de_two',
        paired_sample=FALSE), "Group information table must be a data frame and cannot be data frame with 0 columns and 0 rows.")
    expect_error(as_summarized_experiment(
        two_abundance, two_char_table, group_info="ABC", se_type='de_two',
        paired_sample=FALSE), "Group information table must be a data frame and cannot be data frame with 0 columns and 0 rows.")
})

# Test for group info error input
test_that("as_summarized_experiment can handle group_info error input.", {
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, group_info_twoGroup, se_type=NULL, paired_sample=FALSE),
        "se_type must be one of 'profiling', 'de_two', 'de_multiple', 'ml', or 'corr'.")
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, group_info_twoGroup, se_type='de_two', paired_sample=NULL),
        "paired_sample must be logical value for two group data of differential analysis.")
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, group_info_twoGroup, se_type='de_multiple', paired_sample=NULL),
        "Group information table must contain columns of sample_name, label_name, and group.")
    expect_error(
        as_summarized_experiment(
            multi_abundance, multi_char_table, group_info_multiGroup, se_type='de_multiple', paired_sample=TRUE),
        "paired_sample should be NULL.")
})

# Test for group info error format
test_that("as_summarized_experiment can handle group_info format error.", {
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, group_info_multiGroup, se_type='de_two', paired_sample=FALSE),
        "Group information table must contain columns of sample_name, label_name, group, and pair.")
    expect_error(
        as_summarized_experiment(multi_abundance, multi_char_table, group_info_twoGroup, se_type='de_multiple', paired_sample=NULL),
        "Group information table must contain columns of sample_name, label_name, and group.")
    wrong_name <- group_info_twoGroup %>% dplyr::rename(paired_sample=pair)
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, wrong_name, se_type='de_two', paired_sample=FALSE),
        "Group information table column names must be arranged in order of sample_name, label_name, group, and pair.")
    wrong_name <- group_info_multiGroup %>% dplyr::rename(group_name=group)
    expect_error(
        as_summarized_experiment(multi_abundance, multi_char_table, wrong_name, se_type='de_multiple', paired_sample=NULL),
        "Group information table column names must be arranged in order of sample_name, label_name, and group.")
    expect_error(
        as_summarized_experiment(two_abundance[, 1:10], two_char_table, group_info_twoGroup, se_type='de_two', paired_sample=FALSE),
        "Sample names in group information table 'sample_name' column must be as same as the sample names in lipid abundance data.")
    wrong_group_info <- group_info_twoGroup
    wrong_group_info[2, "group"] <- NA
    expect_error(
        as_summarized_experiment(two_abundance, two_char_table, wrong_group_info, se_type='de_two', paired_sample=FALSE),
        "NA values are not allowed in the 'sample_name', 'label_name', and 'group' columns.")
    wrong_twoGroup <- group_info_twoGroup
    wrong_twoGroup[5:10, "group"] <- "A"
    expect_error(
        as_summarized_experiment(two_abundance, two_char_table, wrong_twoGroup, se_type='de_two', paired_sample=FALSE),
        "The column 'group' in group information table must contain 2 groups.")
    wrong_multiGroup <- group_info_multiGroup
    wrong_multiGroup[14:17, "group"] <- "DHA"
    expect_error(
        as_summarized_experiment(multi_abundance, multi_char_table, wrong_multiGroup, se_type='de_multiple', paired_sample=NULL),
        "The column 'group' in group information table must contain more than 2 groups.")
    wrong_multiGroup <- group_info_multiGroup[group_info_multiGroup$sample_name %in% c("ctrl1", "DHA1", "AA1"), ]
    multi_oneSample <- multi_abundance[, c("feature", "ctrl1", "DHA1",  "AA1")]
    expect_error(
        as_summarized_experiment(multi_oneSample, multi_char_table, wrong_multiGroup, se_type='de_multiple', paired_sample=NULL),
        "In Group information table of multiple groups, each group must contains more than 2 samples.")
    ## pair sample for two group
    expect_error(
        as_summarized_experiment(
            two_abundance, two_char_table, group_info_twoGroup, se_type='de_two', paired_sample=TRUE),
        "Incorrect values in group information table 'pair' column. Each pair
        should be sequentially numbered from 1 to N without any missing, blank,
        or skipped numbers; otherwise, the value must be marked as NA.")
    group_info_twoGroup <- group_info_twoGroup[1:22, ] %>% dplyr::select(-pair) %>%
        dplyr::mutate(pair=1:22)
    expect_error(
        as_summarized_experiment(
            two_abundance[, 1:23], two_char_table, group_info_twoGroup, se_type='de_two', paired_sample=TRUE),
        "Incorrect values in group information table 'pair' column. Each pair
        should be sequentially numbered from 1 to N without any missing, blank,
        or skipped numbers; otherwise, the value must be marked as NA.")
    group_info_twoGroup[2:4, "pair"] <- NA
    expect_error(
        as_summarized_experiment(
            two_abundance[, 1:23], two_char_table, group_info_twoGroup, se_type='de_two', paired_sample=TRUE),
        "Incorrect values in group information table 'pair' column. Each pair
        should be sequentially numbered from 1 to N without any missing, blank,
        or skipped numbers; otherwise, the value must be marked as NA.")
})

# Test for correlation analysis
test_that("as_summarized_experiment can handle wrong correlation group_info.", {
    group_info <- corr_group_info[, 1:2]
    expect_error(
        as_summarized_experiment(
            corr_abundance, corr_char_table, group_info,
            se_type='corr', paired_sample=NULL),
        "Group information table must include columns for sample_name and at least two clinical terms."
    )
    group_info <- corr_group_info[, 1:2] %>% dplyr::mutate(A="A", .after=sample_name)
    expect_error(
        as_summarized_experiment(
            corr_abundance, corr_char_table, group_info,
            se_type='corr', paired_sample=NULL),
        "Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values."
    )
    group_info <- corr_group_info[, 1:3]
    group_info[2:4, 3] <- NA
    expect_error(
        as_summarized_experiment(
            corr_abundance, corr_char_table, group_info,
            se_type='corr', paired_sample=NULL),
        "Except for the first column, sample_name, the group information table must include at least two columns with numeric values and no missing values."
    )
    group_info <- corr_group_info[, 1:3] %>% dplyr::mutate(A="A", .after=sample_name)
    corr_se <- as_summarized_experiment(
        corr_abundance, corr_char_table, corr_group_info, se_type='corr', paired_sample=NULL)
    expect_correct_se(corr_se)
})


# Test for error data frame input
test_that("as_summarized_experiment can handle wrong input.", {
    expect_error(as_summarized_experiment(
        two_char_table, two_abundance, group_info_twoGroup, se_type='de_two',
        paired_sample=FALSE))
    expect_error(as_summarized_experiment(
        group_info_twoGroup, two_char_table, two_abundance, se_type='de_two',
        paired_sample=FALSE))
})

expect_correct_charList <- function(se, result_list) {
    group_info <- as.data.frame(SummarizedExperiment::colData(se))
    expect_true(is.list(result_list), info = "Output should be a list")
    expect_true("deChar_list" %in% names(result_list), info = "'deChar_list' should be in the list")
    expect_true("chain_db_list" %in% names(result_list), info = "'chain_db_list' should be in the list")
    expect_true("common_list" %in% names(result_list), info = "'common_list' should be in the list")
    expect_true("ml_char_list" %in% names(result_list), info = "'ml_char_list' should be in the list")
    expect_true(is.vector(result_list$common_list))
    if (ncol(group_info) > 1) {
        if (all(unique(group_info$group)==c(0,1)) ) {
            expect_true(is.vector(result_list$ml_char_list))
            expect_identical(result_list$deChar_list, "Only provided when differential analysis data is input.")
            expect_identical(result_list$chain_db_list, "Only provided when differential analysis data is input.")
        } else {
            expect_identical(result_list$ml_char_list, "Only provided when machine learning analysis data is input.")
            expect_true(is.vector(result_list$deChar_list))
            expect_true(is.vector(result_list$chain_db_list))
        }
    } else {
        expect_identical(result_list$deChar_list, "Only provided when differential analysis data is input.")
        expect_identical(result_list$chain_db_list, "Only provided when differential analysis data is input.")
        expect_identical(result_list$ml_char_list, "Only provided when machine learning analysis data is input.")
    }
}

expect_correct_seList <- function(se, result_list) {
    group_info <- as.data.frame(SummarizedExperiment::colData(se))
    expect_true(is.list(result_list), info = "Output should be a list")
    expect_true("abundance" %in% names(result_list), info = "'abundance' should be in the list")
    expect_true("lipid_char_table" %in% names(result_list), info = "'lipid_char_table' should be in the list")
    expect_s3_class(result_list$abundance, "data.frame")
    expect_s3_class(result_list$lipid_char_table, "data.frame")
    if (ncol(group_info) > 1) {
        expect_true("group_info" %in% names(result_list), info = "'group_info' should be in the list")
        expect_s3_class(result_list$lipid_char_table, "data.frame")
    }
}

profiling_raw <- as_summarized_experiment(
    two_abundance, two_char_table, group_info=NULL, se_type='profiling', paired_sample=NULL)
profiling_se <- data_process(
    profiling_raw, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
two_raw <- as_summarized_experiment(
    two_abundance, two_char_table, group_info_twoGroup, se_type="de_two", paired_sample=FALSE)
two_se <- data_process(
    two_raw, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
multi_raw <- suppressWarnings(
    as_summarized_experiment(
        multi_abundance, multi_char_table, group_info_multiGroup, se_type="de_multiple", paired_sample=NULL))
multi_se <- data_process(
    multi_raw, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')
ml_raw <- suppressWarnings(
    as_summarized_experiment(
        ml_abundance, ml_char_table, condition_table, se_type='ml', paired_sample=NULL)
    )
ml_se <- data_process(
    ml_raw, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')

corr_raw <- suppressWarnings(
    as_summarized_experiment(
        corr_abundance, corr_char_table, corr_group_info, se_type='corr', paired_sample=NULL)
)
corr_se <- data_process(
    corr_raw, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5, normalization='Percentage', transform='log10')

# Test for basic functionality and output structure
test_that("list_lipid_char can output correct lipid list.", {
    result_list <- list_lipid_char(profiling_se)
    expect_correct_charList(profiling_se, result_list)
    result_list <- list_lipid_char(two_se)
    expect_correct_charList(two_se, result_list)
    result_list <- suppressWarnings( list_lipid_char(multi_se) )
    suppressWarnings( expect_correct_charList(multi_se, result_list) )
    result_list <- list_lipid_char(ml_se)
    expect_correct_charList(ml_se, result_list)
    result_list <- list_lipid_char(corr_se)
    expect_correct_charList(corr_se, result_list)
})

# Test for wrong input
test_that("list_lipid_char can handle wrong input", {
    expect_error(list_lipid_char(data.frame()))
    expect_error(list_lipid_char(abundance_twoGroup))
    expect_error(list_lipid_char(abundance_twoGroup$feature))
})


# Test for basic functionality and output structure
test_that("extract_summarized_experiment can output correct result list", {
    result_list <- extract_summarized_experiment(profiling_se)
    expect_correct_seList(profiling_se, result_list)
    result_list <- extract_summarized_experiment(two_se)
    expect_correct_seList(two_se, result_list)
    expect_error(extract_summarized_experiment(data.frame()))
    expect_error(extract_summarized_experiment(abundance_twoGroup))
    expect_error(extract_summarized_experiment(abundance_twoGroup$feature))
    ## test for have metadata
    deSp_se <- deSp_twoGroup(
        two_se, ref_group='ctrl', test='t-test', significant='padj',
        p_cutoff=0.05, FC_cutoff=2, transform='log10')
    result_list <- extract_summarized_experiment(deSp_se)
    expect_correct_seList(deSp_se, result_list)
})
