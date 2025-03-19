#' @title data_process
#' @description This function processes the abundance data based on user options,
#' including removing features with missing values, imputing missing values, and
#' normalization.
#' @param se A SummarizedExperiment object construct by \code{\link{as_summarized_experiment}}.
#' @param exclude_missing Logical. If exclude_missing=TURE, lipids with
#' missing values will be removed. Default is \code{TRUE}.
#' @param exclude_missing_pct Numeric. Lipids with missing values over a
#' certain percentage (5-100) should be removed. Default is \code{70}.
#' @param replace_na_method Character. The method for NA values replacing.
#' Allowed methods include "QRILC", "SVD", "KNN", "IRMI", "min", "mean", "median",
#' "PPCA", "BPCA", "RandomForest", and "none". If you have already replaced NAs,
#' select 'none'. Default is \code{'min'}.
#' @param replace_na_method_ref Numeric. The value for replacing NA values
#' varies depending on the selected method, and each method applies different number ranges.
#' \enumerate{
#' \item QRILC: 0.1-1
#' \item SVD: 1-10
#' \item KNN: 1-10
#' \item min: 0.1-0.5
#' \item PPCA: 1-10
#' \item BPCA: 1-10
#' }
#' Default is \code{0.5} for replace_na_method='min'.
#' @param normalization Character. Normalization function. Allowed methods
#' include "Percentage", "PQN", "Quantile", "Sum", "Median", and "none".
#' If you have already normalized the abundance values, select 'none'.
#' Default is \code{'Percentage'}.
#' @return Return a SummarizedExperiment object with processed abundance values.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_data <- data_process(de_data_twoGroup, exclude_missing=TRUE,
#'     exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')

data_process <- function(
      se, exclude_missing=TRUE, exclude_missing_pct=70,
      replace_na_method=c('none', 'QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', 'RandomForest'),
      replace_na_method_ref=0.5,
      normalization=c('none', 'Percentage', 'PQN', 'Quantile', 'Sum', 'Median'),
      transform=c('none', 'log10', 'cube', 'square')){
   # check input SE
   .check_inputSE(se, metadata_list=NULL)
   abundance <- .extract_df(se, type = "abundance")
   lipid_char_table <- .extract_df(se, type = "lipid")
   group_info <- .extract_df(se, type="group")
   abundance[-1] <- lapply(abundance[-1], as.numeric)
   # check parameter
   if (!is.logical(exclude_missing)) {
      stop("exclude_missing must be a logical value.")
   } else if (isTRUE(exclude_missing)) {
      if (isFALSE(is.numeric(exclude_missing_pct) && .check_numeric_range(exclude_missing_pct, 5, 100)) ){
         stop("exclude_missing_pct must be a numeric value between 5 and 100.")
      }
   }
   if (isTRUE(any(sapply(abundance[-1], function(x) any(x == 0 | is.na(x)))) &&
              replace_na_method=="none")) {
      stop("Detect 0 values or NAs in abundance data. Please select a data imputation method in replace_na_method other than 'none'.")
   }
   if (is.null(replace_na_method) | isFALSE(replace_na_method %in% c('QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', 'RandomForest')) ) {
      stop("replace_na_method must be one of 'none', 'QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', or 'RandomForest'.")
   } else if (replace_na_method %in% c('QRILC', 'SVD', 'KNN', 'min', 'PPCA', 'BPCA')) {
      if (!is.numeric(replace_na_method_ref)) {
         stop("replace_na_method_ref must be a numeric value.")
      }
   }
   if (is.null(normalization) | isFALSE(normalization %in% c('Percentage', 'PQN', 'Quantile', 'Sum', 'Median', 'none')) ) {
      stop("normalization must be one of 'Percentage', 'PQN', 'Quantile', 'Sum', 'Median', or 'none'.")
   }
   if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'cube', 'square')) ) {
       stop("transfomation must be one of 'log10', 'cube', 'square', or 'none'.")
   }

   # replace 0 to NA
   abundance[abundance==0] <- NA
   # exclude_missing
   if(exclude_missing){
      abundance <- .exclude_missing(abundance, exclude_missing_pct)
   }
   # check if all data has been removed
   if (nrow(abundance)==0){
      stop("All abundance data has been removed, please check and reset the parameter.")
   }
   if(replace_na_method != 'none'){
      abundance <- .missing_impute(
         abundance, replace_na_method=replace_na_method, replace_na_method_ref=replace_na_method_ref)
   }
   if(normalization != 'none'){
      abundance <-  .normalization(abundance, normalization)
   }
   if(transfomation != 'none'){
       processed_abund <-  .transform(abundance, transfomation)
   }


   abundance_mat <- abundance %>% dplyr::arrange(feature) %>%
      tibble::column_to_rownames(var="feature")

   lipid_char_table_trans <- as.data.frame(lipid_char_table[
      (lipid_char_table[[1]] %in% abundance[[1]]),]) %>% dplyr::arrange(feature)
   colnames(lipid_char_table_trans) <- colnames(lipid_char_table)

   transform_SE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(abundance=as.matrix(abundance_mat) ),
      rowData=S4Vectors::DataFrame(lipid_char_table_trans, row.names=lipid_char_table_trans$feature),
      colData=SummarizedExperiment::colData(se),
      metadata=list(processed_abund=processed_abund,
                    transform=transfomation))

   return(transform_SE)
} #function: data_process()

.exclude_missing <- function(abundance, exclude_missing_pct){
   abundance2 <- abundance[,-1]
   missing_pct <- apply(abundance2, 1, function(x){sum(is.na(x))/length(x)})
   maintain_var <- missing_pct*100 < exclude_missing_pct
   abundance <- abundance[maintain_var,]
   return(abundance)
}

.missing_impute <- function(
      abundance, replace_na_method='min', replace_na_method_ref=NULL){
   if (replace_na_method %in% c('QRILC', 'SVD', 'KNN', 'min', 'PPCA', 'BPCA')) {
      if (!is.numeric(replace_na_method_ref)) {
         stop("replace_na_method_ref must be a numeric value.")
      }
   }
   abundance2 <- abundance[-1]
   min_exp <- min(unlist(abundance2)[unlist(abundance2) > 0], na.rm=TRUE)
   if (replace_na_method == 'QRILC'){
      if (.check_numeric_range(replace_na_method_ref, 0.1, 1)) {
         abundance2 <-
            imputeLCMD::impute.QRILC(
               abundance2, tune.sigma=replace_na_method_ref)[[1]]
      } else {
         stop("replace_na_method_ref must be a numeric value between 0.1 and 1.")
      }
   }else if(replace_na_method == 'SVD'){
      if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
         abundance2 <- imputeLCMD::impute.wrapper.SVD(
            abundance2, K=replace_na_method_ref)
      } else {
         stop("replace_na_method_ref must be a numeric value between 1 and 10.")
      }
   }else if(replace_na_method == 'KNN'){
      if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
         abundance2 <- imputeLCMD::impute.wrapper.KNN(
            as.matrix(abundance2), K=replace_na_method_ref)
      } else {
         stop("replace_na_method_ref must be a numeric value between 1 and 10.")
      }
   }else if (replace_na_method == 'IRMI'){
      abundance2_mat <-  as.data.frame(t(VIM::irmi(t(abundance2), imp_var = FALSE)))
      colnames(abundance2_mat) <- colnames(abundance2)
      abundance2 <- abundance2_mat
   }else if(replace_na_method == 'min'){
      if (.check_numeric_range(replace_na_method_ref, 0.1, 0.5)) {
         for (a in seq_len(nrow(abundance2))) {
            abundance2[a,][is.na(abundance2[a,])] <- replace_na_method_ref*min_exp
         }
      } else {
         stop("replace_na_method_ref must be a numeric value between 0.1 and 0.5.")
      }
   }else if(replace_na_method == 'mean'){
      for(a in seq_len(nrow(abundance2))){
         abundance2[a,][is.na(abundance2[a,])] <- mean(
            unlist(abundance2[a,]), na.rm=TRUE)
      }
   }else if(replace_na_method == 'median'){
      for(a in seq_len(nrow(abundance2))){
         abundance2[a,][is.na(abundance2[a,])] <-
            stats::median(unlist(abundance2[a,]), na.rm=TRUE)
      }
   }else if(replace_na_method == 'PPCA'){
      if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
         abundance2 <- pcaMethods::completeObs(
            pcaMethods::pca(abundance2, nPcs=replace_na_method_ref, method="ppca"))
      } else {
         stop("replace_na_method_ref must be a numeric value between 1 and 10.")
      }
   }else if(replace_na_method == 'BPCA'){
      if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
         abundance2 <- pcaMethods::completeObs(
            pcaMethods::pca(abundance2, nPcs=replace_na_method_ref, method="bpca"))
      } else {
         stop("replace_na_method_ref must be a numeric value between 1 and 10.")
      }
   }else if(replace_na_method == 'RandomForest'){
      abundance2 <- missForest::missForest(abundance2)[["ximp"]]
   } else {
      stop("replace_na_method must be one of 'QRILC', 'SVD', 'KNN', 'IRMI', 'min', 'mean', 'median', 'PPCA', 'BPCA', or 'RandomForest'.")
   }
   ## dealing negative values
   numNeg <- sum(abundance2 %>% as.data.frame() %>% dplyr::summarise_all(function(x) sum(x < 0)))
   if (any(numNeg)) {
      warning(paste0(
         numNeg, " negative values were detected after imputation with the ", replace_na_method,
         " method; these values were reverted to half-minute imputation.")
      )
      min_exp <- min(unlist(abundance2)[unlist(abundance2) > 0], na.rm=TRUE)
      abundance2[abundance2<0] <- NA
      for(a in seq_len(nrow(abundance2))){
         abundance2[a,][is.na(abundance2[a,])] <- 0.5*min_exp
      }
   }
   ##
   abundance <- cbind(abundance[1], abundance2)
   return(abundance)
}

.normalization <- function(abundance, normalization='Percentage'){
   if(normalization == 'Percentage'){
      abundance[-1] <-
         purrr::map2(abundance[-1],
                     colSums(abundance[-1], na.rm=TRUE), ~.x/.y*100) %>%
         as.data.frame()
   }else if(normalization == 'PQN'){
      abundance[-1] <- t(Rcpm::pqn(t(abundance[-1])))
   }else if(normalization == 'Quantile'){
      abundance[-1] <- t(
         preprocessCore::normalize.quantiles(t(abundance[-1]), copy=FALSE))
   }else if(normalization == "Sum"){
      abundance[-1] <-
         purrr::map2(abundance[-1],
                     colSums(abundance[-1], na.rm=TRUE), ~.x/.y) %>%
         as.data.frame()
   }else if(normalization == "Median"){
      abundance[-1] <-
         purrr::map2(abundance[-1],
                     colMeans(abundance[-1], na.rm=TRUE), ~.x/.y) %>%
         as.data.frame()
   } else {
      stop("normalization must be one of 'Percentage', 'PQN', 'Quantile', 'Sum', or 'Median'")
   }
   ## check no value is lower than 0
   if (isTRUE(.check_nor_negative(abundance[-1])) ) {
      stop("Detect negative values in normalized abundance data. Please select a different method for normalization.")
   }
   return(abundance)
}

.transform <- function(
        abundance, transform=c('none', 'log10', 'square', 'cube', 'log2')){
    switch(transform,
           log10 = abundance[-1] <- log10(abundance[-1]),
           square = abundance[-1] <- sqrt(abundance[-1]),
           cube = abundance[-1] <- abundance[-1]^(1/3),
           log2 = abundance[-1] <- log2(abundance[-1]))
    return(abundance)
}
