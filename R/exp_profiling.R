#' @title exp_profiling
#' @description This function provides a simple view of sample variability, and
#' compare the amount/expression difference of lipid between samples (i.e.,
#' patients vs. control). The returned tables can be further applied for
#' plotting.
#' @param exp_data_SE A SummarizedExperiment object contains information about
#' various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data  corresponds to specific lipid features, such as class
#' and total length. The first column's name must be "feature" (lipid species),
#' and NAs are allowed for this data. The column data comprises sample names,
#' sample labels, group names, and pair numbers that represent 'the pair' for
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of the rowData of \bold{exp_data_SE},
#' such as total length.
#' @return Return a list with 3 SummarizedExperiment objects.
#' \enumerate{
#' \item tot.num.lip: number of expressed lipids and total amount of lipid in
#' each sample.
#' \item dens.lip: the underlying probability distribution of the lipid
#' expression in each sample.
#' \item exp.compo.by.lipid: the expression level of each sample within each
#' group (e.g., PE, PC) of selected characteristics (e.g., class).
#' }
#' @export
#' @examples
#' library(dplyr)
#' library(SummarizedExperiment)
#' data("profiling_data")
#' exp_data_SE <- profiling_data
#' lipid_char_table <- as.data.frame(rowData(exp_data_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' exp_profiling(exp_data_SE, char_var[1])
exp_profiling <- function(exp_data_SE, char_var=NULL){
  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c("feature",
                          SummarizedExperiment::colData(exp_data_SE)[[1]])
  lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))

   exp_trans_data <- exp_data %>%
    tidyr::gather(sample_name, value, -feature)
  num.lip <- exp_trans_data %>%
    dplyr::mutate(is.expr=!is.na(value)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(expr.count=sum(is.expr)) %>%
    dplyr::arrange(sample_name)
  num.lip$sample_name <- factor(
    num.lip$sample_name, levels=unique(num.lip$sample_name)[order(
      num.lip$expr.count, decreasing=TRUE)])
  tot.lip <- exp_trans_data %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(lipid_amount=sum(value, na.rm=TRUE))
  tot.lip$sample_name <- factor(
    tot.lip$sample_name, levels=unique(tot.lip$sample_name)[order(
      tot.lip$lipid_amount, decreasing=TRUE)])
  dens.lip <- stats::na.omit(exp_trans_data)
  if(!is.null(lipid_char_table)){
    lipid_char_table <- as.data.frame(lipid_char_table)
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
        comma_count <- max(stringr::str_count(FA_lipid_char_table[, col],
                                              ','), na.rm=TRUE)
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
          select_name <- paste0("_", i)
          FA_lipid_char_table <-FA_lipid_char_table[-intersect(grep(
            select_name, FA_lipid_char_table[, "lipid.category"]),
            which(is.na(FA_lipid_char_table$lipid.category.value))),]
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
    FA_col <- grep("FA_", colnames(lipid_char_table), value=TRUE)
    if(length(FA_col) > 0){
      max_comma <- 0
      for(i in seq_len(length(FA_col))){
        col <- FA_col[i]
        comma_count <- max(stringr::str_count(lipid_char_table[, col],
                                              ','), na.rm=TRUE)
        if(comma_count > 0){
          lipid_char_table <- lipid_char_table %>%
            tidyr::separate(col, c(col, paste0(col, "_", seq_len(comma_count))))
        }
        if(comma_count > max_comma){max_comma <- comma_count}
      }
      lipid_char_table <- lipid_char_table %>%
        tidyr::gather(lipid.category, lipid.category.value, -feature)
      if(max_comma > 0){
        for (i in seq_len(max_comma)) {
          select_name <- paste0("_",i)
          lipid_char_table <-lipid_char_table[-intersect(grep(
            select_name, lipid_char_table[,"lipid.category"]),
            which(is.na(lipid_char_table$lipid.category.value))),]
        }
        for(i in seq_len(length(FA_col))){
          col <- FA_col[i]
          lipid_char_table[grep(col,
                                lipid_char_table[, "lipid.category"]),
                           "lipid.category"] <- col
        }
      }
    }else{
      lipid_char_table <- lipid_char_table %>%
        tidyr::gather(lipid.category, lipid.category.value, -feature)
    }
    exp.compo.by.lipid <- merge(exp_data %>%
                                  tidyr::gather(sample_name, value, -feature),
                                lipid_char_table, by="feature") %>%
      dplyr::filter(lipid.category == char_var &
                      !is.na(lipid.category.value)) %>%
      dplyr::group_by(sample_name, lipid.category, lipid.category.value) %>%
      dplyr::summarise(value=sum(value, na.rm=TRUE)) %>%
      plotly::ungroup() %>%
      dplyr::group_by(sample_name) %>%
      dplyr::mutate(weight=100/sum(value)) %>%
      dplyr::mutate(value=value*weight)%>%
      plotly::ungroup()
  }else{
    exp.compo.by.lipid = NULL
  }

  tot.num.lip <- merge(x = num.lip, y = tot.lip, by = "sample_name", all = TRUE)

  tot.num.lip_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(tot.num.lip=tot.num.lip))

  dens.lip_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(dens.lip=dens.lip))

  exp.compo.by.lipid_SE <- SummarizedExperiment::SummarizedExperiment(
    assays=list(exp.compo.by.lipid=exp.compo.by.lipid))

  return(list(tot.num.lip=tot.num.lip_SE, dens.lip=dens.lip_SE,
              exp.compo.by.lipid=exp.compo.by.lipid_SE))
}
