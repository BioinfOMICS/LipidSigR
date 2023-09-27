#' @title corr_heatmap
#' @description The correlation heatmap illustrates the correlation between
#' lipid species or samples and depicts the patterns in each group.
#' The correlation can be calculated by Pearson or Spearman.
#' \enumerate{
#' \item Compute and cluster the correlation coefficient depending on
#' user-defined method and distance.
#' \item Print two correlation heatmaps to illustrate the correlation between
#' lipid species or samples.
#' }
#' @param exp_data_SE A SummarizedExperiment object contains information about
#' various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data  corresponds to specific lipid features, such as class
#' and total length. The first column's name must be "feature" (lipid species),
#' and NAs are allowed for this data. The column data comprises sample names,
#' sample labels, group names, and pair numbers that represent 'the pair' for
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param corr_method A character string indicating which correlation
#' coefficient is to be computed. One of "pearson" or
#' "spearman". (default: "pearson")
#' @param distfun A character string of the distance measure indicating which
#' correlation coefficient (or covariance) is to be computed. Allowed methods
#' include "pearson", "spearman", "kendall", euclidean", "maximum",
#' "manhattan", "canberra", "binary", "minkowski". (default: "maximum")
#' @param hclustfun A character string of the agglomeration method to be used.
#' This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). (default: "average")
#' @param type It should be 'sample' or 'lipid.' 'sample' outputs the
#' correlation results of samples, and 'lipid' outputs output the correlation
#' results of lipid species.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' data("profiling_data")
#' exp_data_SE <- profiling_data
#' exp_transform_SE <- data_process(exp_data_SE, exclude_var_missing=TRUE,
#'     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
#'     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
#'     data_transform=TRUE, trans_type='log', centering=FALSE, scaling=FALSE)
#' corr_heatmap(exp_transform_SE, corr_method="pearson", distfun="maximum",
#'     hclustfun="average",type='sample')

corr_heatmap <- function(exp_data_SE, corr_method="pearson", distfun="maximum",
                         hclustfun='average', type='sample'){

  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature",
                          SummarizedExperiment::colData(exp_data_SE)[[1]])

  cor.test.p <- function(x,method){
    FUN <- function(x, y,method) stats::cor.test(x=x, y=y,
                                                 method=method,
                                                 na.action=
                                                   "na.exclude")[["p.value"]]
    z <- outer(colnames(x), colnames(x),
               Vectorize(function(i, j) FUN(x[, i], x[, j], method=method)))
    dimnames(z) <- list(colnames(x), colnames(x))
    z
  }
  if(type == 'sample'){
    exp_data <- subset(exp_data, select=-feature)
    corr_coef <- stats::cor(exp_data, method=corr_method,
                            use="pairwise.complete.obs")
    corr_p <- cor.test.p(x=exp_data, method=corr_method)
    data <- corr_coef
  }else if(type == 'lipid'){
    exp_data <- exp_data %>%
      tidyr::gather(-feature, key='sample_name', value='value') %>%
      tidyr::spread(key='feature', value='value')
    exp_data <- subset(exp_data, select=-sample_name)
    corr_coef <- stats::cor(exp_data, method=corr_method,
                            use="pairwise.complete.obs")
    corr_p <- cor.test.p(x=exp_data, method=corr_method)
    data <- corr_coef
  }

  heatmap_data <- SummarizedExperiment::SummarizedExperiment(
    assays=list(corr_coef=corr_coef, corr_p=corr_p))

  return(heatmap_data)
}
