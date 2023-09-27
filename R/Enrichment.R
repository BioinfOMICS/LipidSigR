#' @title Enrichment
#' @description
#' Determine whether significant lipid species are enriched in
#' user-selected characteristics.
#' @param DE_species_sig A SummarizedExperiment object comprises the significant
#' results of differential expression analysis, including fold change, p-value,
#' adjusted p-value. The output of \code{\link{DE_species}}.
#' @param exp_transform_SE A SummarizedExperiment object contains information
#' about various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' not allowed. The row data corresponds to specific lipid features, such as
#' class and total length. The first column's name must be "feature" (lipid
#' species), and NAs are allowed for this data. The column data comprises
#' sample names, sample labels, group names, and pair numbers that represent
#' 'the pair' for conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of the rowData of \bold{exp_data_SE},
#' such as total length.
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @importFrom magrittr %>%
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
#' Enrichment(DE_result_sig, exp_transform_SE, char_var=char_var[1],
#'     sig_pvalue=0.05)
Enrichment <- function(DE_species_sig, exp_transform_SE, char_var,
                       sig_pvalue=0.05){

  DE_species_table_sig <- as.data.frame(
    SummarizedExperiment::assay(DE_species_sig))
  colnames(DE_species_table_sig)[1] <- 'feature'

  lipid_char_table <- as.data.frame(
    SummarizedExperiment::rowData(exp_transform_SE))

  if(ncol(dplyr::select(lipid_char_table,
                        tidyselect::starts_with("FA_"))) == 0){
    warning("(OPTIONAL) lipid_char_table does not contain column
            names starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>%
      dplyr::select(feature, tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_", colnames(FA_lipid_char_table), value=TRUE)
    max_comma <- 0
    for(i in seq_len(length(FA_col))){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[, col],
                                            ','),na.rm=TRUE)
      if(comma_count > 0){
        FA_lipid_char_table <- FA_lipid_char_table %>%
          tidyr::separate(col,
                          c(col,
                            paste0(col, "_", seq_len(comma_count))),
                          ",",
                          convert=TRUE)
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>%
      tidyr::gather(lipid.category, lipid.category.value, -feature)
    if(max_comma > 0){
      for (i in seq_len(max_comma)) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <-FA_lipid_char_table[-intersect(
          grep(select_name,
               FA_lipid_char_table[, "lipid.category"]),
          which(is.na(FA_lipid_char_table$lipid.category.value))), ]
      }
    }
    if(methods::is(FA_lipid_char_table$lipid.category.value, 'character')){
      stop("In the 'FA_' related analyses, the values are positive integer
           or zero and separated by comma. i.e., 10,12,11")
    }else if(sum(
      stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)) !=
      round(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value)))) != 0 |
      min(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value))) < 0){
      stop("In the 'FA_' related analyses, the values are positive integer
           or zero and separated by comma. i.e., 10,12,11")
    }
  }
  lipo.fisher <- function(sig.mat, lipo.category){
    sig.class <- sig.mat %>%
      dplyr::group_by(characteristic)%>%
      dplyr::summarise(sig.count=dplyr::n())%>%
      dplyr::left_join(lipo.category, by='characteristic')%>%
      dplyr::mutate(p.value=NA, significant='NO')
    n.sig <- sig.class %>% dplyr::summarise(sum=sum(sig.count)) %>% .$sum
    n.lipid <- lipo.category %>% dplyr::summarise(sum=sum(total.count)) %>%
      .$sum
    for (i in seq_len(nrow(sig.class))){
      sig.class$p.value[i] <- stats::fisher.test(
        matrix(c(sig.class$sig.count[i], sig.class$total.count[i],
                 n.sig-sig.class$sig.count[i],
                 n.lipid-sig.class$total.count[i]),
               nrow=2), alternative='greater')$p.value
    }
    sig.class$m.log.p <- -log10(sig.class$p.value)
    sig.class$significant[sig.class$p.value < sig_pvalue] <- 'YES'
    return(sig.class)
  }
  lipo.category <- lipid_char_table %>%
    dplyr::select(feature, 'characteristic'=tidyselect::all_of(char_var)) %>%
    dplyr::group_by(characteristic) %>%
    dplyr::summarise(total.count=dplyr::n())
  char.tab <- lipid_char_table %>%
    dplyr::select(feature, 'characteristic'=tidyselect::all_of(char_var))
  DE_species_table_sig <- DE_species_table_sig %>%
    dplyr::left_join(char.tab, by='feature') %>%
    dplyr::filter(!is.na(characteristic))
  sig_lipid.all <- DE_species_table_sig %>%
    dplyr::mutate(condition='UP & DOWN')
  sig_lipid.up <- DE_species_table_sig %>%
    dplyr::filter(log2FC > 0) %>%
    dplyr::mutate(condition='UP')
  sig_lipid.down <- DE_species_table_sig %>%
    dplyr::filter(log2FC < 0) %>%
    dplyr::mutate(condition='DOWN')
  sig_class.all <- lipo.fisher(sig_lipid.all, lipo.category) %>%
    dplyr::mutate(condition='UP & DOWN')
  if(nrow(sig_lipid.up) > 0){
    sig_class.up <- lipo.fisher(sig_lipid.up, lipo.category) %>%
      dplyr::mutate(condition='UP')
  }else{
    sig_class.up <- NULL
  }
  if(nrow(sig_lipid.down) > 0){
    sig_class.down <- lipo.fisher(sig_lipid.down, lipo.category) %>%
      dplyr::mutate(condition='DOWN')
  }else{
    sig_class.down <- NULL
  }
  sig_class <- dplyr::bind_rows(sig_class.all, sig_class.up, sig_class.down)
  sig_class$characteristic <- factor(sig_class$characteristic,
                                     sort(unique(sig_class$characteristic)))
  final.tab <- sig_class %>%
    dplyr::filter(condition != 'UP & DOWN') %>%
    dplyr::mutate(mlogP=ifelse(condition == 'DOWN', -(m.log.p), m.log.p),
        significance=ifelse(significant == 'YES' & condition == 'UP', 'UP',
            ifelse(significant == 'YES' & condition == 'DOWN',
                   'DOWN', 'non-significant')))
  final.tab$m.log.p <- round(final.tab$m.log.p, 3)

  final.tab_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(final.tab=as.data.frame(final.tab)))

  return(final.tab_SE)
}
