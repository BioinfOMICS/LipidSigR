#' @title DE_char
#' @description
#' The \code{\link{DE_char}} function is for lipid characteristics analysis,
#' which is used to search for deferentially expressed lipid characters based
#' on the user's choice. In this function, two-way ANOVA will be used to
#' evaluate the interaction between lipid characteristics (class, total length,
#' etc.) and group information (control v.s experimental group) on the lipid
#' expression values or percentages (according to selected characteristics). To
#' obtain the overview of differentially expressed lipids by user-selected
#' characteristics, an expression table from \code{\link{Species2Char}} and a
#' lipid characteristic table are needed. Note that the sample can be paired or
#' not, depending on \bold{paired=TRUE/FALSE}.
#' @param exp_data_SE A SummarizedExperiment object contains information
#' about various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data corresponds to specific lipid features, such as class
#' and total length. The first  column's name must be "feature" (lipid
#' species), and NAs are allowed for this data. The column data comprises
#' sample names, sample labels, group names, and  pair numbers that represent
#' 'the pair' for conducting t-tests or Wilcoxon  tests. NAs are allowed.
#' \emph{NOTE: The output of Species2Char}.
#' @param data_transform Logical. If data_transform = TRUE, transformed
#'     exp_data by log10.
#' @param paired Logical. If \bold{paired = TRUE}, data are paired samples.
#'     (default: FALSE)
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @param sig_FC Numeric. Significance of the fold-change. (default: 2)
#' @param insert_ref_group A character string. The name of 'ctrl' after name
#' conversion.
#' @param ref_group A character string. The name of 'exp' after name conversion.
#' @return Return a list with 5 SummarizedExperiment objects.
#' \enumerate{
#' \item char_exp_data: the expression of the user-selected lipid
#' characteristics by samples.
#' \item char_table_all: statistical analysis for differentially expressed
#' characters. For example, the result of two-way ANOVA will tell you whether
#' the interaction effect between totallength and different groups is present.
#' The post hoc test (t-test) will then calculate significance for each length
#' of FA and produce a p-value.
#' \item combined_table: the value calculated by the weighted average of the
#' expression of the continuous lipid characteristics.
#' \item combine_result_table: statistics of t.test for control and experiment
#' groups.
#' \item bar_table: data that can be further applied for plotting bar plot.
#' \item char_box: data that can be further applied for plotting box plot.
#' \item char_table_all_sig: the significant results of statistical analysis
#' for differentially expressed characters.
#' }
#' @export
#' @examples
#' library(SummarizedExperiment)
#' data("DE_data")
#' exp_data_SE <- DE_data
#' lipid_char_table <- as.data.frame(rowData(exp_data_SE))
#' char_var <- colnames(lipid_char_table)[-1]
#' Spe2Char_result <- Species2Char(exp_data_SE, char_var = char_var[4])
#' exp_transform_non_log <- data_process(Spe2Char_result,
#'     exclude_var_missing=TRUE, missing_pct_limit=50, replace_zero=TRUE,
#'     zero2what='min', xmin=0.5, replace_NA=TRUE, NA2what='min', ymin=0.5,
#'     pct_transform=TRUE, data_transform=FALSE, trans_type='log',
#'     centering=FALSE, scaling=FALSE)
#' DE_char(exp_transform_non_log, data_transform = TRUE, paired = FALSE,
#'     sig_pvalue = 0.05, sig_FC = 2, insert_ref_group=NULL, ref_group=NULL)
DE_char <- function(exp_data_SE, data_transform=TRUE, paired=FALSE,
    sig_pvalue=0.05, sig_FC=2, insert_ref_group=NULL, ref_group=NULL){

  exp_data <- cbind(SummarizedExperiment::rowData(exp_data_SE)[[1]],
                    as.data.frame(SummarizedExperiment::assay(exp_data_SE)))
  colnames(exp_data) <- c(
    colnames(SummarizedExperiment::rowData(exp_data_SE))[[1]],
    SummarizedExperiment::colData(exp_data_SE)[[1]] )
  group_info <- as.data.frame(SummarizedExperiment::colData(exp_data_SE))
  rownames(group_info) <- NULL

  char_exp_data <- exp_data
  var <- char_exp_data[[1]]
  char_name <- colnames(char_exp_data)[1]
  if(data_transform == TRUE){
    Char_Transform_table <- cbind(char_exp_data[1], log10(char_exp_data[-1]))
  }else{
    Char_Transform_table <- char_exp_data
  }
  colnames(Char_Transform_table)[1] <- 'feature'
  if(length(unique(Char_Transform_table[[1]])) > 1){
    two_way_anova <- Char_Transform_table %>%
      tidyr::gather(-feature, key='sample_name', value='value') %>%
      dplyr::left_join(group_info[c(1, 3)], by='sample_name') %>%
      dplyr::mutate(feature=as.factor(feature))
    two_way_anova <- tryCatch({stats::aov(value ~ group+feature+group:feature,
                                          data=two_way_anova)},
                              error=function(e){NULL})
    if(is.null(two_way_anova)){
      anova_pvalue <- NA
    }else{
      if(nrow(summary(two_way_anova)[[1]]) != 4){
        anova_pvalue <- NA
      }else{
        anova_pvalue <- summary(two_way_anova)[[1]][[3, "Pr(>F)"]]
      }
    }
  }else{
    anova_pvalue <- NA
  }
  group1 <- group_info %>% dplyr::filter(group == 'ctrl') %>%
    dplyr::arrange(pair) %>% .$sample_name
  group2 <- group_info %>% dplyr::filter(group == 'exp') %>%
    dplyr::arrange(pair) %>% .$sample_name
  var_name <- character(length(var))
  mean_ctrl <- numeric(length(var))
  sd_ctrl <- numeric(length(var))
  mean_exp <- numeric(length(var))
  sd_exp <- numeric(length(var))
  FC <- numeric(length(var))
  pvalue <- numeric(length(var))
  for(a in seq_len(length(var))){
    var_name[a] <- var[a]
    mean_ctrl[a] <- mean(unlist(char_exp_data[a, group1]), na.rm=TRUE)
    mean_exp[a] <- mean(unlist(char_exp_data[a, group2]), na.rm=TRUE)
    sd_ctrl[a] <- stats::sd(unlist(char_exp_data[a, group1]), na.rm=TRUE)
    sd_exp[a] <- stats::sd(unlist(char_exp_data[a, group2]), na.rm=TRUE)
    FC[a] <- mean(unlist(char_exp_data[a, group2]), na.rm=TRUE)/
      mean(unlist(char_exp_data[a, group1]), na.rm=TRUE)
    if(paired == TRUE){
      pvalue[a] <- tryCatch(
        stats::t.test(unlist(Char_Transform_table[a, group1]),
                      unlist(Char_Transform_table[a, group2]),
                      paired=TRUE, var.equal=TRUE)$p.value,
        error=function(e){NA}
      )
    }else{
      pvalue[a] <- tryCatch(
        stats::t.test(unlist(Char_Transform_table[a, group1]),
                      unlist(Char_Transform_table[a, group2]),
                      paired=FALSE, var.equal=TRUE)$p.value,
        error=function(e){NA}
      )
    }
  }
  Result_table <- data.frame(var_name=var_name, method='two-way anova',
                             anova_pvalue=anova_pvalue,
                             post_hoc_test='t.test', mean_ctrl=mean_ctrl,
                             sd_ctrl=sd_ctrl, mean_exp=mean_exp,
                             sd_exp=sd_exp, FC=FC, log2FC=log2(FC),
                             post_hoc_pvalue=pvalue)
  Result_table <- Result_table %>%
    dplyr::mutate(sig=ifelse(post_hoc_pvalue < sig_pvalue &
                               abs(log2FC) > log2(sig_FC), 'yes', 'no'))
  Result_table <- Result_table %>%
    dplyr::mutate(sig=ifelse(is.na(sig), 'no',sig ))
  Result_table[is.na(Result_table)] <- NA
  colnames(Result_table)[1] <- char_name
  if(sum(is.na(as.numeric(var))) == 0){
    char_exp_data[[1]] <- as.numeric(char_exp_data[[1]])
    Combined_char_data <- purrr::map2(char_exp_data[1],
                                      char_exp_data[-1],
                                      ~sum(.x*.y,
                                           na.rm = TRUE)/sum(.y,
                                                             na.rm=TRUE)) %>%
      as.data.frame() %>%
      dplyr::mutate(Combine_name=stringr::str_c(char_name, '_index')) %>%
      dplyr::select(Combine_name, dplyr::everything())
    colnames(Combined_char_data) <- colnames(char_exp_data)
    mean_ctrl <- mean(unlist(Combined_char_data[, group1]), na.rm=TRUE)
    mean_exp <- mean(unlist(Combined_char_data[, group2]), na.rm=TRUE)
    sd_ctrl <- stats::sd(unlist(Combined_char_data[, group1]), na.rm=TRUE)
    sd_exp <- stats::sd(unlist(Combined_char_data[, group2]), na.rm=TRUE)
    FC <- mean_exp/mean_ctrl
    if(data_transform == TRUE){
      Combined_char_transform_table <- cbind(Combined_char_data[1],
                                             log10(Combined_char_data[-1]))
    }else{
      Combined_char_transform_table <- Combined_char_data
    }
    if(paired == TRUE){
      pvalue <- tryCatch(
        stats::t.test(unlist(Combined_char_transform_table[, group1]),
                      unlist(Combined_char_transform_table[, group2]),
                      paired=TRUE,
                      var.equal=TRUE)$p.value,
        error=function(e){NA}
      )
    }else{
      pvalue <- tryCatch(
        stats::t.test(unlist(Combined_char_transform_table[, group1]),
                      unlist(Combined_char_transform_table[, group2]),
                      paired=FALSE,
                      var.equa=TRUE)$p.value,
        error=function(e){NA}
      )
    }
    Combine_char_result_table <- data.frame(var_name=Combined_char_data[1,1],
                                            method='t.test',
                                            mean_ctrl=mean_ctrl,
                                            sd_ctrl=sd_ctrl,
                                            mean_exp=mean_exp,
                                            sd_exp=sd_exp,
                                            FC=FC, log2FC=log2(FC),
                                            p_value=pvalue)
    Combine_char_result_table <- Combine_char_result_table %>%
      dplyr::mutate(sig=ifelse(p_value < sig_pvalue &
                                 abs(log2FC) > log2(sig_FC), 'yes', 'no'))
    Combine_char_result_table <- Combine_char_result_table %>%
      dplyr::mutate(sig=ifelse(is.na(sig), 'no', sig ))
    Combine_char_result_table[is.na(Combine_char_result_table)] <- NA
    colnames(Combine_char_result_table)[1] <- colnames(Combined_char_data)[1]
  }else{
    Combined_char_data <- data.frame()
    Combine_char_result_table <- data.frame()
  }
  CHAR <- colnames(Result_table[1])
  CTRL.RES <- Result_table %>%
    dplyr::select(1, sig, mean_ctrl, sd_ctrl) %>%
    dplyr::mutate(Group='Ctrl')
  colnames(CTRL.RES) <- c('Category', 'Significant', 'Mean', 'SD', 'Group')

  EXP.RES <- Result_table %>%
    dplyr::select(1, sig, mean_exp, sd_exp) %>%
    dplyr::mutate(Group='Exp')
  colnames(EXP.RES) <- c('Category', 'Significant', 'Mean', 'SD', 'Group')
  barTab <- data.table::rbindlist(l=list(CTRL.RES, EXP.RES),
                                  use.names=TRUE, fill=TRUE)
  barTab <- barTab %>%
    dplyr::group_by(Category) %>%
    dplyr::mutate(max_error_bar=max(Mean+SD)) %>% dplyr::ungroup()
  barTab$post_hoc_pvalue <- NA
  for(i in seq_len(nrow(barTab))){
    barTab$post_hoc_pvalue[i] <- Result_table$post_hoc_pvalue[
      which(barTab$Category[i] == Result_table[, 1])]
  }
  if(!is.null(insert_ref_group) & !is.null(ref_group)){
    exp_raw_name <- ref_group[-which(insert_ref_group == ref_group)]
    barTab$Group[which(barTab$Group == 'Ctrl')] <-  insert_ref_group
    barTab$Group[which(barTab$Group == 'Exp')] <-  exp_raw_name
    barTab$Group <- factor(barTab$Group,
                           levels=c(insert_ref_group, exp_raw_name))
  }
  if(sum(is.na(as.numeric(barTab$Category))) == 0){
    barTab$Category <- as.factor(as.numeric(barTab$Category))
  }
  if(nrow(Combined_char_data) > 0 & nrow(Combine_char_result_table) > 0){
    boxTab <- Combined_char_data %>%
      tibble::column_to_rownames(var=CHAR) %>%
      t() %>% as.data.frame() %>%
      merge(group_info, by.x=0, by.y='sample_name')
    colnames(boxTab)[2] <- 'Category'
    if(!is.null(insert_ref_group) & !is.null(ref_group)){
      exp_raw_name <- ref_group[-which(insert_ref_group == ref_group)]
      boxTab$group[which(boxTab$group == 'ctrl')] <-  insert_ref_group
      boxTab$group[which(boxTab$group == 'exp')] <-  exp_raw_name
    }

    boxTab_SE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(boxTab=as.data.frame(boxTab)))
  }else{
    boxTab_SE <- NULL
  }

  ## significant table
  char_table_all_sig <- Result_table %>% dplyr::filter(sig=="yes")

  result_list <- c("char_exp_data", "Result_table", "Combined_char_data",
      "Combine_char_result_table", "barTab", "char_table_all_sig")
  for (i in seq(result_list)){
    save_data <- SummarizedExperiment::SummarizedExperiment(
      assays=list(as.data.frame(get0(result_list[i]))))
    assign(paste0(result_list[i], "_SE"), save_data)
  }

  return(list(char_exp_data=char_exp_data_SE,
              char_table_all=Result_table_SE,
              combined_table=Combined_char_data_SE,
              combine_result_table=Combine_char_result_table_SE,
              bar_table=barTab_SE,
              char_box=boxTab_SE,
              char_table_all_sig=char_table_all_sig_SE))
}
