#' @title dim_process
#' @description This function transforms the data into a format suitable for
#' further use in dimension reduction methods, such as PCA, t-SNE, UMAP,
#' and PLSDA.
#' @param exp_transform_SE A SummarizedExperiment object contains information
#' about various features, such as molecules, lipid class, etc., and their
#' expression levels for each sample. The assay is a matrix representing the
#' expression of lipid features across all samples. Missing values (NAs)  are
#' allowed. The row data  corresponds to specific lipid features, such as class
#' and total length. The first column's name must be "feature" (lipid species),
#' and NAs are allowed for this data. The column data comprises sample names,
#' sample labels, group names, and pair numbers that represent 'the pair' for
#' conducting t-tests or Wilcoxon tests. NAs are allowed.
#' @param sig_feature A character string for identifying significant variable in
#' differential analysis results. It is only required when \bold{type}
#' is "PLSDA."
#' @param type A character of \bold{"PCA"}, \bold{"tsne"}, \bold{"UMAP"},
#' and \bold{"PLSDA"}. The dimension reduction method to which the transformed
#' data will be applied.
#' @param insert_ref_group A character string. The name of 'ctrl' after name
#' conversion.
#' @param ref_group A character string. The name of 'exp' after name conversion.
#' @return Return a SummarizedExperiment object.
#' @export
#' @examples
#' data("DE_data")
#' exp_data_SE <- DE_data
#' exp_transform_SE <- data_process(exp_data_SE, exclude_var_missing=TRUE,
#'     missing_pct_limit=50, replace_zero=TRUE, zero2what='min', xmin=0.5,
#'     replace_NA=TRUE, NA2what='min', ymin=0.5,  pct_transform=TRUE,
#'     data_transform=FALSE, trans_type='log', centering=FALSE, scaling=FALSE)
#' dim_process(exp_transform_SE, sig_feature=NULL,
#'     type='PCA', insert_ref_group=NULL, ref_group=NULL)
dim_process <- function(exp_transform_SE, sig_feature=NULL, type='PLSDA',
                        insert_ref_group=NULL, ref_group=NULL){

  exp_transform_table <- cbind(
    SummarizedExperiment::rowData(exp_transform_SE)[[1]],
    as.data.frame(SummarizedExperiment::assay(exp_transform_SE)))
  colnames(exp_transform_table) <- c(
    colnames(SummarizedExperiment::rowData(exp_transform_SE))[[1]],
    SummarizedExperiment::colData(exp_transform_SE)[[1]] )

  group_info <- as.data.frame(SummarizedExperiment::colData(exp_transform_SE))
  rownames(group_info) <- NULL

  if(!is.null(sig_feature)){
    exp_transform_table <- exp_transform_table %>%
      dplyr::filter(eval(
        parse(text=colnames(exp_transform_table)[1])) %in% sig_feature)
  }
  num <- apply(exp_transform_table[-1], 1, FUN=function(x){length(unique(x))})
  exp_transform_table <- exp_transform_table[(num!=1),]
  exp_transform_table <- exp_transform_table[!is.infinite(rowSums(
    exp_transform_table[-1], na.rm=TRUE)),]
  if(type == 'PLSDA'){
    if(!is.null(group_info)){
      group_info <- group_info %>%
        dplyr::filter(sample_name %in% colnames(exp_transform_table)[-1])
      if(!is.null(insert_ref_group) & !is.null(ref_group) &
         !is.null(group_info)){
        Y <- ifelse(group_info$group == insert_ref_group, 1, 2)
      }else{
        Y <- ifelse(group_info$group == 'ctrl', 1, 2)
      }
      names(Y) <- group_info$sample_name
      rownames(exp_transform_table) <- exp_transform_table$feature
      dim_process_table <- exp_transform_table[
        c(colnames(SummarizedExperiment::rowData(exp_transform_SE))[[1]],
                                                 group_info$sample_name)] %>%
        dplyr::select(-1) %>%
        t()
    }else{
      stop('group_info can not be NULL')
    }
  }else{
    dim_process_table <- exp_transform_table[-1] %>% t() %>% as.data.frame()
    colnames(dim_process_table) <- exp_transform_table[[1]]
    dim_process_table <- stats::na.omit(dim_process_table)
    Y <- NULL
  }

  dim_data <- SummarizedExperiment::SummarizedExperiment(
    assays=list(dim_process_table=dim_process_table), rowData=Y)

  return(dim_data)
}
