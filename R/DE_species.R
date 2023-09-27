#' @title DE_species
#' @description Compute differentially expressed analysis of
#' two groups (independent) based on the group Information to find
#' significant lipid species. To find differentially expressed lipid species,
#' a statistical method (t-test or Wilcoxon test) and the cut-offs for
#' significant lipid species need to be chosen, and the p-value will then be
#' adjusted by Benjamini-Hochberg or other methods.
#' @param exp_data_SE A SummarizedExperiment object contains information about 
#' various features, such as molecules, lipid class, etc., and their 
#' expression levels for each sample. The assay is a matrix representing the 
#' expression of lipid features across all samples. Missing values (NAs)  are 
#' allowed. The row data  corresponds to specific lipid features, such as class 
#' and total length. The first column's name must be "feature" (lipid species), 
#' and NAs are allowed for this data. The column data comprises sample names, 
#' sample labels, group names, and pair numbers that represent 'the pair' for 
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param data_transform Logical. If data_transform = TRUE,
#' transform exp_data by log10.
#' @param paired Logical. If paired = TRUE,
#' data are paired samples. (default: FALSE)
#' @param test A character string indicating which method to be used for
#' comparing means. Allowed method include \bold{"t.test"} and
#' \bold{"wilcox.test"}. (default: "t.test")
#' @param adjust_p_method Correction method, a character string. One of
#' \bold{"holm"}, \bold{"hochberg"}, \bold{"hommel"}, \bold{"bonferroni"},
#' \bold{"BH"}, \bold{"BY"}, \bold{"fdr"}, \bold{"none"}, can be
#' abbreviated. (default: "BH")
#' @param sig_stat A character string indicating which p-value is to be used
#' for the statistically significant. One of \bold{"p.adj"} or
#' \bold{"p"}. (default: "p.adj")
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @param sig_FC Numeric. Significance of the fold-change. (default: 2)
#' @return Return a list with 2 SummarizedExperiment objects.
#' \enumerate{
#' \item exp_data_stat: an SE object comprises the results of
#' differential expression analysis, including fold change, p-value,
#' adjusted p-value.
#' \item exp_data_stat_sig: an SE object comprises the significant
#' results of differential expression analysis, including fold change,
#' p-value, adjusted p-value.
#' }
#' @export
#' @examples
#' data("DE_data")
#' exp_data_SE <- DE_data
#' exp_transform_non_log <- data_process(exp_data_SE, exclude_var_missing=TRUE,
#'     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5, 
#'     replace_NA=TRUE, NA2what='min', ymin=0.5, pct_transform=TRUE,
#'     data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
#' DE_species(exp_data_SE=exp_transform_non_log, data_transform=FALSE, 
#'     paired=FALSE, test='wilcoxon test', adjust_p_method='BH',
#'     sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2)
DE_species <- function(exp_data_SE, data_transform=TRUE, paired=FALSE, 
                       test='t.test', adjust_p_method='BH',
                       sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2){
  
  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]], 
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature", 
                          SummarizedExperiment::colData(exp_data_SE)[[1]])
  
  group_info <- as.data.frame(SummarizedExperiment::colData(exp_data_SE))
  rownames(group_info) <- NULL
  
  exp_data_ga <- exp_data %>%
    tidyr::gather(sample_name, value, -1) %>%
    dplyr::left_join(group_info, by='sample_name')
  if(paired == TRUE){
    exp_data_ga <- exp_data_ga %>%
      dplyr::arrange(group, pair, feature)
  }
  if(data_transform == TRUE){
    exp_data_trans_ga <- exp_data_ga
    exp_data_trans_ga[,3] <- log10(exp_data_trans_ga[, 3])
    exp_data_tab <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
                    mean_exp=mean(unlist(exp), na.rm=TRUE),
                    method='t-test',
                    FC=mean_exp/mean_ctrl,
                    log2FC=log2(FC))
    exp_data_p <- exp_data_trans_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature)
    if(test == 't.test'){
      exp_data_tab <- exp_data_tab %>% dplyr::mutate(method='t-test')
      exp_data_p <- exp_data_p %>%
        dplyr::mutate(p_value=tryCatch(
          stats::t.test(unlist(ctrl), unlist(exp),
                        paired=paired, var.equal=TRUE)$p.value,
          error=function(e){NA}))
    }else{
      exp_data_tab <- exp_data_tab %>% dplyr::mutate(method='wilcoxon test')
      exp_data_p <- exp_data_p %>%
        dplyr::mutate(p_value=tryCatch(
          stats::wilcox.test(unlist(ctrl), unlist(exp), paired=paired)$p.value,
          error=function(e){NA}))
    }
    exp_data_p <- exp_data_p %>% dplyr::mutate(m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value, 
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)
    exp_data_stat <- exp_data_tab[-2:-3] %>%
      dplyr::right_join(exp_data_p[-2:-3], by='feature')
  }else if(data_transform == FALSE){
    exp_data_p <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
                    mean_exp=mean(unlist(exp), na.rm=TRUE),
                    FC=mean_exp/mean_ctrl,
                    log2FC=log2(FC))
    if(test == 't.test'){
      exp_data_p <- exp_data_p %>%
        dplyr::mutate(method='t-test',
                      p_value=tryCatch(
                        stats::t.test(unlist(ctrl), unlist(exp),
                                      paired=paired, var.equal=TRUE)$p.value,
                        error=function(e){NA}))
    }else{
      exp_data_p <- exp_data_p %>%
        dplyr::mutate(method='wilcoxon test',
                      p_value=tryCatch(
                        stats::wilcox.test(unlist(ctrl), unlist(exp),
                                           paired=paired)$p.value,
                        error=function(e){NA}))
    }
    exp_data_p <- exp_data_p %>% dplyr::mutate(m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value,
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)
    exp_data_stat <- exp_data_p[-2:-3]
  }
  exp_data_stat <- exp_data_stat %>%
    dplyr::mutate(sig_p=ifelse(p_value < sig_pvalue &
                                 abs(log2FC) > log2(sig_FC), 'yes', 'no'),
                  sig_p_adj=ifelse(p_adj < sig_pvalue &
                                     abs(log2FC) > log2(sig_FC), 'yes', 'no'))
  exp_data_stat$mean_ctrl <- round(exp_data_stat$mean_ctrl, 5)
  exp_data_stat$mean_exp <- round(exp_data_stat$mean_exp, 5)
  exp_data_stat$FC <- round(exp_data_stat$FC, 3)
  exp_data_stat$log2FC <- round(exp_data_stat$log2FC, 3)
  exp_data_stat$m_log10_p_value <- round(exp_data_stat$m_log10_p_value, 3)
  exp_data_stat$m_log10_p_adj <- round(exp_data_stat$m_log10_p_adj, 3)
  
  ## significant 
  if(sum(!is.na(exp_data_stat[,'p_value'])) == 0){
    sig.diff.exp <- exp_data_stat %>%
      dplyr::filter(abs(log2FC) >= log2(sig_FC)) %>%
      dplyr::distinct(feature, .keep_all=TRUE)
    warning('Only two samples could not be tested, significant lipids were
            instead defined by log2 fold change.')
  }else{
    if(sig_stat == 'p'){
      sig.diff.exp <- exp_data_stat %>%
        dplyr::filter(sig_p == 'yes') %>%
        dplyr::distinct(feature, .keep_all=TRUE)
    }else if(sig_stat == 'p.adj'){
      sig.diff.exp <- exp_data_stat %>%
        dplyr::filter(sig_p_adj == 'yes') %>%
        dplyr::distinct(feature, .keep_all=TRUE)
    }
  }
  
  
  exp_data_stat_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(exp_data_stat=as.data.frame(exp_data_stat)))
  
  exp_data_stat_sig <- SummarizedExperiment::SummarizedExperiment(
    assays=list(exp_data_stat_sig=as.data.frame(sig.diff.exp)))
  
  
  return(list(exp_data_stat=exp_data_stat_SE,
              exp_data_stat_sig=exp_data_stat_sig))
}