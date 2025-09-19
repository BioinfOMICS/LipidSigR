library(testthat)
# library(SummarizedExperiment)
# library(tidyverse)
# library(LipidSigR)

data("abundance_multiGroup")

return.objects <- function(results){
    # Check class of return object
    expect_s3_class(results, 'data.frame')
    # Check NA
    expect_true(all(!is.na(results[, c(1:4)])))
    # Check numeric
    expect_true(all(sapply(results[, c(11:13)], is.numeric)))
    # Check number of columns
    expect_equal(ncol(results), 72)
}


test_that("Test under normal conditions", {
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_multiGroup$feature)
    recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
    goslin_annotation <- parse_lipid %>%
        dplyr::filter(Original.Name %in% recognized_lipid)
    results <- lipid_annotation(goslin_annotation)
    return.objects(results)
})

test_that("Test when none of the lipids can be recognized", {
    abundance_multiGroup$feature <- paste0(abundance_multiGroup$feature, '!@#$')
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_multiGroup$feature[1:20])
    expect_error(lipid_annotation(parse_lipid), info='None of the lipid names are parseable by Goslin. Please check if your lipid dialects conform to the formats supported by Goslin.')
})

test_that("Test with only 'Sn.Position' level lipids", {
    abundance_multiGroup %<>% dplyr::filter(stringr::str_detect(feature, '_'))
    abundance_multiGroup$feature <- stringr::str_replace_all(
        abundance_multiGroup$feature, pattern='_', replacement='\\/')
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=sample(abundance_multiGroup$feature, 20))
    recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
    goslin_annotation <- parse_lipid %>%
        dplyr::filter(Original.Name %in% recognized_lipid)
    results <- lipid_annotation(goslin_annotation)
    return.objects(results)
})

test_that("Test with only long-chain base (LCB) lipids", {
    abundance_multiGroup %<>% .[1:2,]
    abundance_multiGroup$feature <- c('SPB 18:1;O2', 'SPB 18:0;O2')
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_multiGroup$feature)
    recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
    goslin_annotation <- parse_lipid %>%
        dplyr::filter(Original.Name %in% recognized_lipid)
    results <- lipid_annotation(goslin_annotation)
    return.objects(results)
})

test_that("Test lipids that cannot be mapped to any resource ID", {
    abundance_multiGroup %<>% .[1:3,]
    abundance_multiGroup$feature <- c('BMP 15:1_16:4', 'BMP 15:2_22:6', 'BMP 15:4_22:5')
    parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_multiGroup$feature)
    recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
    goslin_annotation <- parse_lipid %>%
        dplyr::filter(Original.Name %in% recognized_lipid)
    results <- lipid_annotation(goslin_annotation)
    return.objects(results)
})



