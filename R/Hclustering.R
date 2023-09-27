#' @title Hclustering
#' @description Hierarchical clustering of lipid species or characteristics
#' derived from two groups.
#' @param exp_data_SE A SummarizedExperiment object contains information
#' about various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' not allowed. The row data corresponds to specific lipid features, such as
#' class and total length. The first column's name must be "feature" (lipid
#' species), and NAs are allowed for this data. The column data comprises
#' sample names, sample labels, group names, and pair numbers that represent
#' 'the pair' for conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param DE_result_sig A SummarizedExperiment object comprises the significant
#' results of differential expression analysis, including fold change, p-value,
#' adjusted p-value. The output of \code{\link{DE_species}}.
#' @param group_info A data frame comprises the name of the sample, the label
#' of the sample, the group name of the sample, and the pair number represents
#' 'the pair' for the t-test/Wilcoxon test. NAs are allowed.
#' @param type It should be 'all' or 'sig.' 'all' for output the results of all
#' lipid species or characteristics, and 'sig' for output the results of
#' significant lipid species or characteristics.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' library(dplyr)
#' library(SummarizedExperiment)
#' data("DE_data")
#' exp_data_SE <- DE_data
#' exp_transform_SE <- data_process(exp_data_SE, exclude_var_missing=TRUE,
#'     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
#'     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
#'     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
#' exp_transform_non_log <- data_process(exp_data_SE, exclude_var_missing=TRUE,
#'     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
#'     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
#'     data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
#' DE_species_result <- DE_species(exp_transform_non_log, data_transform=TRUE,
#'     paired=FALSE, test='t.test', adjust_p_method='BH', sig_stat='p.adj',
#'     sig_pvalue=0.05, sig_FC=2)
#' lipid_char_table <- as.data.frame(rowData(exp_transform_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' DE_result_sig <- DE_species_result$exp_data_stat_sig
#' Hclustering(exp_data_SE=exp_transform_SE,
#'     DE_result_sig=DE_result_sig, type='all')
Hclustering <- function(exp_data_SE, DE_result_sig, group_info, type='all'){

  exp_data <- cbind(
    SummarizedExperiment::rowData(exp_data_SE)[[1]],
    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature",
    SummarizedExperiment::colData(exp_data_SE)[[1]])
  rownames(exp_data) <- NULL

  group_info <- as.data.frame(SummarizedExperiment::colData(exp_data_SE))
  rownames(group_info) <- NULL

  DE_result_table <- as.data.frame(SummarizedExperiment::assay(DE_result_sig))
  colnames(DE_result_table)[1] <- 'feature'

  exp.mat <- exp_data %>%
    dplyr::select(feature, group_info$sample_name)
  if(type=='sig'){
    exp.mat <- exp.mat %>%
      dplyr::filter(feature %in% DE_result_table$feature)
  }
  exp.mat <- exp.mat %>%
    tibble::column_to_rownames(var='feature') %>%
    as.matrix()
  colnames(exp.mat) <- group_info$label_name
  exp.mat <- sweep(exp.mat, 1, rowMeans(exp.mat, na.rm=TRUE))
  exp.mat <- sweep(exp.mat, 1, apply(exp.mat, 1, sd, na.rm=TRUE), "/")
  if(sum(is.na(exp.mat)) > 0){
    exp.mat <- exp.mat[-which(is.na(exp.mat[, 2])), ]
  }

  exp.mat_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(exp.mat=exp.mat))

  return(exp.mat_SE)
}
