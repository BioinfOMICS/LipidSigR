#' @title Sig_lipid_feature
#' @description \code{\link{Sig_lipid_feature}} shows the significant
#' lipid species based on different lipid characteristics and visualizes the
#' difference between control and experimental groups by applying
#' log2 Fold Change.
#' @param DE_result_sig A SummarizedExperiment object comprises the significant
#' results of differential expression analysis, including fold change, p-value,
#' adjusted p-value. The output of \code{\link{DE_species}}.
#' @param exp_transform_SE A SummarizedExperiment object contains information
#' about various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' not allowed. The row data corresponds to specific lipid features, such as
#' class and total length. The first column's name must be "feature" (lipid
#' species), and NAs are allowed for this data. The column data comprises sample
#' names, sample labels, group names, and pair numbers that represent
#' 'the pair' for conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of the rowData of \bold{exp_data_SE},
#' such as total length.
#' @param sig_FC Numeric. Significance of the fold-change. (default: 2)
#' @return Return a list with 2 SummarizedExperiment object.
#' @export
#' @examples
#' library(magrittr)
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
#' Sig_lipid_feature(DE_result_sig, exp_transform_SE, char_var[1], sig_FC=2)
Sig_lipid_feature <- function(DE_result_sig, exp_transform_SE,
                              char_var, sig_FC=2){

  DE_species_table_sig <- as.data.frame(
    SummarizedExperiment::assay(DE_result_sig))
  colnames(DE_species_table_sig)[1] <- 'feature'

  lipid_char_table <- as.data.frame(
    SummarizedExperiment::rowData(exp_transform_SE))

  if(ncol(dplyr::select(lipid_char_table,
                        tidyselect::starts_with("FA_"))) == 0){
    warning("(OPTIONAL) lipid_char_table does not contain column names
            starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>%
      dplyr::select(feature, tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_", colnames(FA_lipid_char_table), value=TRUE)
    max_comma <- 0
    for(i in seq_len(length(FA_col))){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[, col], ','),
                         na.rm=TRUE)
      if(comma_count > 0){
        FA_lipid_char_table <- FA_lipid_char_table %>%
          tidyr::separate(col, c(col, paste0(col, "_", seq_len(comma_count))),
                          ",", convert=TRUE)
      }
      if(comma_count > max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>%
      tidyr::gather(lipid.category, lipid.category.value, -feature)
    if(max_comma > 0){
      for (i in seq_len(max_comma)) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <-FA_lipid_char_table[-intersect(
          grep(select_name, FA_lipid_char_table[, "lipid.category"]),
          which(is.na(FA_lipid_char_table$lipid.category.value))), ]
      }
    }
    if(methods::is(FA_lipid_char_table$lipid.category.value, 'character')){
      stop("In the 'FA_' related analyses, the values are positive integer or
           zero and separated by comma. i.e., 10,12,11")
    }else if(sum(
      stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)) !=
      round(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value)))) != 0 |
      min(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value))) < 0){
      stop("In the 'FA_' related analyses, the values are positive integer or
           zero and separated by comma. i.e., 10,12,11")
    }
  }
  char <- lipid_char_table %>%
    dplyr::select(feature, tidyselect::all_of(char_var))
  colnames(char)[2] <- 'characteristic'
  plot.tab <- DE_species_table_sig %>%
    dplyr::left_join(char, by='feature') %>%
    dplyr::filter(!is.na(characteristic)) %>%
    dplyr::arrange(characteristic)
  plot.tab$characteristic <- factor(plot.tab$characteristic,
                                    sort(unique(plot.tab$characteristic)))
  plot.tab <- plot.tab %>% dplyr::group_by(characteristic)
  sig_class.log2FC <- plot.tab %>%
    dplyr::summarise(log2FC.mean=mean(log2FC), log2FC.sd=sd(log2FC),
                     log2FC.direction=log2FC.mean / abs(log2FC.mean),
                     freqs=dplyr::n(),
                     significant='NO')%>%
    dplyr::mutate(log2FC.meansd=log2FC.mean+log2FC.sd*log2FC.direction,
                  significant=replace(x=significant,
                                      list=abs(log2FC.mean) > log2(sig_FC),
                                      values='YES'))
  sig.dotchart <- plot.tab %>% dplyr::mutate(log2FC.mean=mean(log2FC))
  sig.dotchart$characteristic <- forcats::fct_reorder(
    sig.dotchart$characteristic, sig.dotchart$log2FC.mean, min)

  sig.class_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(sig_class.log2FC=as.data.frame(sig_class.log2FC)))

  sig.dotchart_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(sig.dotchart=as.data.frame(sig.dotchart)))

  return(list(sig.class_SE=sig.class_SE,
              sig.dotchart_SE=sig.dotchart_SE))
}
