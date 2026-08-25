#' @title data_process
#' @description This function processes the abundance data based on user options,
#' including removing features with missing values, imputing missing values, and
#' normalization.
#' @param se A SummarizedExperiment object construct by \code{\link{as_summarized_experiment}}.
#' @param exclude_missing Logical. If exclude_missing=TURE, lipids with
#' missing values will be removed. Default is \code{TRUE}.
#' @param exclude_missing_method Character. The method for filtration type.
#' Allowed methods include "global", "group".
#' @param exclude_missing_pct Numeric. Lipids with missing values over a
#' certain percentage (5-100) should be removed. Default is \code{70}.
#' @param replace_na_method Character. The method for NA values replacing.
#' Allowed methods include "min", "BPCA", "GMS", "KNN", "impSeq", "Impseqrob",
#' "LLS", "mean", "median", "MinProb", "PPCA", "PI", "QRILC", "RandomForest",
#' "SVD", and "none". If you have already replaced NAs, select "none".
#' Default is \code{'min'}.
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
#' include "Percentage", "CubicSplines", "CyclicLoess", "EigenMS", "linear",
#' "Mean", "Median", "nonLinear", "PQN", "Quantile", "RLR", "Sum", "TMM", and
#' "none". If you have already normalized the abundance values, select 'none'.
#' Default is \code{'Percentage'}.
#' @param transform Character. Transformation function. Allowed methods
#' include "log10", "log2", "cube", "square", and "none".
#' If you have already transformed the abundance values, select 'none'.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object with processed abundance values.
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_data <- data_process(de_data_twoGroup, exclude_missing=TRUE,
#'     exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')

data_process <- function(
      se, exclude_missing=TRUE, exclude_missing_method='global',
      exclude_missing_pct=70, replace_na_method='min', replace_na_method_ref=0.5,
      normalization='Percentage', transform='log10'){
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
              replace_na_method == "none")) {
      stop("Detect 0 values or NAs in abundance data. Please select a data imputation method in replace_na_method other than 'none'.")
   }
   match.arg(replace_na_method,
             c("min", "BPCA", "GMS", "KNN", "impSeq", "Impseqrob", "LLS",
               "mean", "median", "MinProb", "PPCA", "PI", "QRILC",
               "RandomForest", "SVD", 'none'))
   if (replace_na_method %in% c('QRILC', 'SVD', 'KNN', 'min', 'PPCA', 'BPCA')) {
       if (!is.numeric(replace_na_method_ref)) {
           stop("replace_na_method_ref must be a numeric value.")
       }
   }
   match.arg(normalization,
             c('Percentage', 'CubicSplines', 'CyclicLoess', 'EigenMS', 'linear',
               'Mean', 'Median', 'nonLinear', 'PQN', 'Quantile',
               'RLR', 'Sum', 'TMM', 'none'))
   match.arg(transform, c('log10', 'log2', 'cube', 'square', 'none'))

   # replace 0 to NA
   abundance[abundance == 0] <- NA
   # exclude_missing
   if(exclude_missing){
      abundance <- .exclude_missing(abundance, group_info, exclude_missing_method, exclude_missing_pct)
   }
   # check if all data has been removed
   if (nrow(abundance)==0){
      stop("All abundance data has been removed, please check and reset the parameter.")
   }
   if(replace_na_method != 'none' & any(is.na(abundance))){
      abundance <- .missing_impute(
         abundance, replace_na_method=replace_na_method, replace_na_method_ref=replace_na_method_ref)
   }
   if(normalization != 'none'){
      abundance <-  .normalization(abundance, group_info, normalization)
   }
   if(transform != 'none'){
       processed_abund <-  .transform(abundance, transform)
   }else{
       processed_abund <-  abundance
   }

   abundance_mat <- abundance %>% dplyr::arrange(feature) %>%
      tibble::column_to_rownames(var="feature") %>%
       dplyr::select(dplyr::all_of(group_info$sample_name))
   lipid_char_table_trans <- as.data.frame(lipid_char_table[
      (lipid_char_table[[1]] %in% abundance[[1]]),]) %>% dplyr::arrange(feature)
   colnames(lipid_char_table_trans) <- colnames(lipid_char_table)
   colData <- se %>% SummarizedExperiment::colData() %>% as.data.frame() %>%
       tibble::remove_rownames()

   transform_SE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(abundance=as.matrix(abundance_mat) ),
      rowData=S4Vectors::DataFrame(lipid_char_table_trans, row.names=lipid_char_table_trans$feature),
      colData=colData,
      metadata=list(processed_abund=processed_abund,
                    transform=transform))

   return(transform_SE)
} #function: data_process()


.exclude_missing <- function(abundance, groupInfo=NULL, method='global', param=80){
    filtrationFeature <- abundance %>%
        tidyr::gather(key='sample_name', value='abundance', -feature) %>%
        dplyr::group_by(feature)
    if(method == 'group'){
        filtrationFeature %<>% merge(groupInfo) %>% dplyr::group_by(feature, group)
    }
    filtrationFeature %<>%
        dplyr::summarise(missing_pct=mean(is.na(abundance))*100) %>%
        dplyr::filter(missing_pct < param) %>%
        dplyr::pull(feature)
    abundance %<>% dplyr::filter(feature %in% filtrationFeature)
    if (nrow(abundance) == 0){
        cli::cli_abort(c("All abundance data has been removed.",
                         "i"="Please check and reset the parameter."))
    }
    return(abundance)
}

.missing_impute <- function(
      abundance, replace_na_method='min', replace_na_method_ref=NULL){
   abundance2 <- abundance[-1]
   min_exp <- min(unlist(abundance2)[unlist(abundance2) > 0], na.rm=TRUE)
   if(replace_na_method == 'BPCA'){
       if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
           abundance2 %<>% t() %>%
               pcaMethods::pca(nPcs=replace_na_method_ref, method="bpca") %>%
               pcaMethods::completeObs() %>% t()
       } else {
           stop("replace_na_method_ref must be a numeric value between 1 and 10.")
       }
   }else if(replace_na_method == 'GMS'){
       abundance2 %<>% .gmsImpute()
   }else if(replace_na_method == 'KNN'){
       if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
           abundance2 %<>% as.matrix() %>%
               imputeLCMD::impute.wrapper.KNN(., K=replace_na_method_ref)
       } else {
           stop("replace_na_method_ref must be a numeric value between 1 and 10.")
       }
   }else if(replace_na_method == 'impSeq'){
       abundance2 %<>% rrcovNA::impSeq()
   }else if(replace_na_method == 'Impseqrob'){
       abundance2 %<>% as.matrix() %>% rrcovNA::impSeqRob() %>%
           magrittr::extract2('x') %>% as.data.frame()
   }else if(replace_na_method == 'LLS'){
       if(replace_na_method_ref >= ncol(abundance2)){
           replace_na_method_ref <- ncol(abundance2) - 1
           message("The replace_na_method_ref must be less than the sample size, so it is replaced by the sample size minus one.")
       }
       abundance2 %<>% t() %>%
           pcaMethods::llsImpute(k=replace_na_method_ref, allVariables=TRUE) %>%
           pcaMethods::completeObs() %>% t() %>% as.data.frame()
   }else if(replace_na_method == 'mean'){
       rowMean <- rowMeans(abundance2, na.rm=TRUE)
       naPositions <- is.na(abundance2)
       abundance2[naPositions] <- rowMean[row(abundance2)[naPositions]]
   }else if(replace_na_method == 'median'){
       rowMedian <- apply(abundance2, 1, median, na.rm=TRUE)
       naPositions <- is.na(abundance2)
       abundance2[naPositions] <- rowMedian[row(abundance2)[naPositions]]
   }else if(replace_na_method == 'min'){
       if (.check_numeric_range(replace_na_method_ref, 0.1, 0.5)) {
           for (a in seq_len(nrow(abundance2))) {
               abundance2[a,][is.na(abundance2[a,])] <- replace_na_method_ref*min_exp
           }
       } else {
           stop("replace_na_method_ref must be a numeric value between 0.1 and 0.5.")
       }
   }else if(replace_na_method == 'MinProb'){
       abundance2 %<>% imputeLCMD::impute.MinProb()
   }else if(replace_na_method == 'PPCA'){
       if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
           success <- FALSE
           attempt <- 0
           while (!success && attempt < 20) {
               attempt <- attempt + 1
               try({
                   abundance2 <- pcaMethods::completeObs(
                       pcaMethods::pca(
                           abundance2, nPcs=replace_na_method_ref, method="ppca"))
                   success <- TRUE
               }, silent=TRUE)
           }
           if (!success) {
               stop("ppca failed after 20 attempts.")
           }
       } else {
           stop("replace_na_method_ref must be a numeric value between 1 and 10.")
       }
   }else if(replace_na_method == 'PI'){
       for(i in 1:nrow(abundance2)){
           temp <- as.numeric(abundance2[i, ])
           if(sum(is.na(temp)) > 0){
               temp.sd <- 0.3 * sd(temp[!is.na(temp)], na.rm = TRUE)
               temp.mean <- mean(temp[!is.na(temp)], na.rm = TRUE) - 1.8 * sd(temp[!is.na(temp)], na.rm = TRUE)
               n.missing <- sum(is.na(temp))
               temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
               abundance2[i, ] <- temp
           }
       }
   }else if (replace_na_method == 'QRILC'){
       if (.check_numeric_range(replace_na_method_ref, 0.1, 1)) {
           abundance2 %<>%
               imputeLCMD::impute.QRILC(tune.sigma=replace_na_method_ref) %>%
               magrittr::extract2(1)
       } else {
           stop("replace_na_method_ref must be a numeric value between 0.1 and 1.")
       }
   }else if(replace_na_method == 'RandomForest'){
       abundance2 <- missForest::missForest(abundance2)[["ximp"]]
   }else if(replace_na_method == 'SVD'){
       if (.check_numeric_range(replace_na_method_ref, 1, 10)) {
           abundance2 %<>% imputeLCMD::impute.wrapper.SVD(K=replace_na_method_ref)
       } else {
           stop("replace_na_method_ref must be a numeric value between 1 and 10.")
       }
   # }else if (replace_na_method == 'IRMI'){
   #    abundance2_mat <-  as.data.frame(t(VIM::irmi(t(abundance2), imp_var = FALSE)))
   #    colnames(abundance2_mat) <- colnames(abundance2)
   #    abundance2 <- abundance2_mat
   }
   ## dealing negative values
   numNeg <- sum(abundance2 %>% as.data.frame() %>% dplyr::summarise_all(function(x) sum(x < 0)))
   if (any(numNeg)) {
      warning(numNeg, " negative values were detected after imputation with the ",
          replace_na_method,
         " method; these values were reverted to half-minute imputation.")
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

.gmsImpute <-function(abundance, alpha=1, nfolds=10){

    row.missing <- rowSums(is.na(abundance))
    if(sum(row.missing==0) < 5){
        return(abundance)
        cli::cli_abort("Need at least 5 compounds without missing peaks.")
    }

    if(is.null(rownames(abundance))){
        rownames(abundance) <- paste('Var',seq_len(dim(abundance)[1]),sep='')
    }
    if(is.null(colnames(abundance))){
        colnames(abundance) <- paste('sample',seq_len(dim(abundance)[2]),sep='')
    }
    miss <- abundance[(rowSums(!is.na(abundance)) >= 3 &
                           rowSums(is.na(abundance)) != 0),]
    miss_large3 <- abundance[(rowSums(!is.na(abundance))<3),]
    nomiss <- abundance[rowSums(is.na(abundance)) == 0,]
    imputed <- miss
    for(i in seq_len(dim(miss)[1])){
        cvfit <- glmnet::cv.glmnet(
            t(nomiss[,!is.na(miss[i,])]), t(miss[i,][!is.na(miss[i,])]),
            alpha=alpha, nfolds=nfolds, family='gaussian')
        imputed[i,][is.na(miss[i,])] <- stats::predict(
            cvfit, newx=t(nomiss[,is.na(miss[i,])]),s="lambda.min")
    }
    imputed.second <- imputed
    for(i in seq_len(dim(miss)[1])){
        imputed.all <- rbind(nomiss,imputed[-i,])
        cvfit <- glmnet::cv.glmnet(
            t(imputed.all[,!is.na(miss[i,])]), t(miss[i,][!is.na(miss[i,])]),
            alpha=alpha, nfolds=nfolds, family='gaussian')
        imputed.second[i,][is.na(miss[i,])] <- stats::predict(
            cvfit, newx=t(imputed.all[,is.na(miss[i,])]), s="lambda.min")
    }
    imputed.final <- rbind(rbind(nomiss,imputed.second),miss_large3)
    imputed.final <- imputed.final[rownames(abundance),]
    return(imputed.final)
}

.normalization <- function(abundance, groupInfo=NULL, normalization='Percentage'){
    if(normalization == 'CubicSplines'){
        abundance[-1] %<>% log() %>% .CubicSplinesNorm() %>% exp()
    }else if(normalization == 'CyclicLoess'){
        abundance[-1] %<>% log() %>% limma::normalizeCyclicLoess() %>% exp()
    }else if(normalization == 'EigenMS'){
        abundance[-1] %<>% log()
        abundance %<>% .EigenMSNorm(., groupInfo) %>% tibble::remove_rownames()
        abundance[-1] %<>% exp()
    }else if(normalization == 'linear'){
        abundance[-1] %<>% .linearNorm()
    }else if(normalization == 'Mean'){
        abundance[-1] %<>% sweep(2, apply(., 2, mean, na.rm=TRUE), "/")
    }else if(normalization == 'Median'){
        abundance[-1] %<>% sweep(2, apply(., 2, stats::median, na.rm=TRUE), "/")
    }else if(normalization == 'nonLinear'){
        abundance[-1] %<>% log() %>% .nonLinearNorm() %>% exp()
    }else if(normalization == 'Percentage'){
        abundance[-1] %<>% sweep(2, colSums(., na.rm=TRUE), "/") %>% `*`(100)
    }else if(normalization == 'PQN'){
        abundance[-1] %<>% t() %>% Rcpm::pqn() %>% t()
    }else if(normalization == 'Quantile'){
        abundance[-1] %<>% as.matrix() %>%
            preprocessCore::normalize.quantiles(copy=FALSE)
    }else if(normalization == 'RLR'){
        abundance[-1] %<>% log() %>% .RLRNorm() %>% exp()
    }else if(normalization == 'Sum'){
        abundance[-1] %<>% sweep(2, colSums(., na.rm=TRUE), "/")
    }else if(normalization == 'TMM'){
        abundance[-1] %<>% edgeR::DGEList(counts=.) %>%
            edgeR::calcNormFactors() %>% edgeR::cpm()
    }

   ## check no value is lower than 0
   if (isTRUE(.check_nor_negative(abundance[-1])) ) {
      stop("Detect negative values in normalized abundance data. Please select a different method for normalization.")
   }
   return(abundance)
}

.CubicSplinesNorm <- function(abundance) {
    splineAbundance <- abundance %>%
        affy::normalize.qspline(samples=0.1, target=apply(abundance, 1, mean))
    rownames(splineAbundance) <- rownames(abundance)
    colnames(splineAbundance) <- colnames(abundance)
    return(splineAbundance)
}

.EigenMSNorm <- function(abundance, groupInfo) {
    rawNames <- colnames(abundance)
    prot.info <- abundance %>% dplyr::select(feature) %>%
        tibble::rowid_to_column('rowid')
    ints_eig1 <- ProteoMM::eig_norm1(
        m=abundance[-1], treatment=as.factor(groupInfo$group),
        prot.info=prot.info)
    if(ints_eig1$h.c != 0){
        ints_norm1 <- ProteoMM::eig_norm2(rv=ints_eig1)
        abundance <- ints_norm1$normalized %>% dplyr::select(-rowid)
    }
    colnames(abundance) <- rawNames
    return(abundance)
}

.linearNorm <- function(abundance) {
    scalingfactors <- mean(apply(abundance, 1, median))/apply(abundance, 2, mean)
    return(t(t(abundance)*scalingfactors))
}

.RLRNorm <- function(abundance) {
    abundance <- as.matrix(abundance)
    sampleMedian <- matrixStats::rowMedians(as.matrix(abundance), na.rm=TRUE)
    .calculateRLMForCol <- function(colIndex, sampleMedian, abundance) {
        lrFit <- MASS::rlm(as.matrix(abundance[, colIndex])~sampleMedian,
            na.action=stats::na.exclude)
        coeffs <- lrFit$coefficients
        coefIntercept <- coeffs[1]
        coefSlope <- coeffs[2]
        globalFittedRLRCol <- (abundance[, colIndex]-coefIntercept)/coefSlope
        globalFittedRLRCol
    }
    RLR_data <- vapply(
        seq_len(ncol(abundance)),
        .calculateRLMForCol,
        rep(0, nrow(abundance)),
        sampleMedian=sampleMedian,
        abundance=abundance)
    colnames(RLR_data) <- colnames(abundance)
    return(RLR_data)
}

.nonLinearNorm <- function(abundance) {
    average.intensity <- apply(abundance, 2, mean)
    median.number <- round(ncol(abundance)/2 + 0.1)
    ordering <- order(average.intensity)
    median.sample.number <- ordering[median.number]
    median.sample <- abundance[, median.sample.number]
    liwong.abundance <- vector()
    .normalize.invariantset <- function(data, ref, prd.td=c(0.003, 0.007)) {
        np <- length(data)
        r.ref <- rank(ref)
        r.array <- rank(data)
        prd.td.adj <- prd.td*10
        i.set <- rep(TRUE, np)
        ns <- sum(i.set)
        ns.old <- ns+50+1
        while ( (ns.old-ns) > 50 ) {
            air <- (r.ref[i.set] + r.array[i.set]) / (2*ns)
            prd <- abs(r.ref[i.set] - r.array[i.set]) / ns
            threshold <- (prd.td.adj[2]-prd.td[1]) * air + prd.td.adj[1]
            i.set[i.set] <- (prd < threshold)
            ns.old <- ns
            ns <- sum(i.set)
            if (prd.td.adj[1] > prd.td[1])
                prd.td.adj <- prd.td.adj * 0.9
        }
        env <- environment()
        tryCatch({
            n.curve <- stats::smooth.spline(ref[i.set], data[i.set])},
            error=function(e) {
                cli::cli_warn("The error is caused by {.fn smooth.spline}, so
                              replace the {.arg cv} parameter with
                              {.val TRUE}.")
                tryCatch({
                    n.curve <- stats::smooth.spline(ref[i.set],
                        data[i.set], cv=TRUE)
                    assign('n.curve', n.curve, envir=env)
                },
                    error=function(e) {
                        stop()
                    }
                )
            }
        )
        return(list(n.curve=n.curve, i.set=i.set))
    }
    for(i in seq_len(ncol(abundance))){
        liwong.model <- .normalize.invariantset(data=abundance[, i],
            ref=median.sample,
            prd.td=c(0.003, 0.007))
        liwong.sample <- stats::predict(liwong.model$n.curve$fit, abundance[,i])
        liwong.abundance <- cbind(liwong.abundance, liwong.sample$y)
    }
    return(liwong.abundance)
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
