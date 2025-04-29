#' @title ml_model
#' @description This function constructs the machine learning model, the output object
#' can be used as input for plotting and further analyses.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character list. Lipid characteristics selected from the ml_char list
#' returned by \code{\link{list_lipid_char}}. Select 'none' to exclude all lipid characteristics.
#' @param ranking_method Character. The ranking method to be computed.
#' Allowed methods include 'p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM',
#' 'Lasso', 'Ridge', 'ElasticNet'. Default is \code{'Random_forest'}.
#' @param ml_method Character. The machine learning method to be computed.
#' Allowed methods include 'Random_forest', 'SVM', 'Lasso', 'Ridge',
#' 'ElasticNet', 'xgboost'. Default is \code{'Random_forest'}.
#' @param split_prop Numeric. The proportion of data to be retained for modeling/analysis.
#' The range is 0.1 to 0.5. Default is \code{0.3}.
#' @param nfold Numeric. The number of fold that the original dataset is randomly
#' partitioned into equal-sized subsamples. Must be a positive interger.
#' Default is \code{10}.
#' @param alpha Numeric. The alpha value between 0 and 1 when choosing
#' \code{"ElasticNet"} as ml_method. 0 is for Ridge and 1 is for Lasso.
#' If "ElasticNet" is not selected as the ml_method, set the value to NULL.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @return Return a SummarizedExperiment object containing analysis results.
#' @export
#' @examples
#' data("ml_sub")
#' processed_se <- data_process(
#'     ml_sub, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' char_list <- list_lipid_char(processed_se)
#' ml_se <- ml_model(
#'     processed_se, char=c("class","Total.DB"), transform='log10',
#'     ranking_method='Random_forest', ml_method='Random_forest', split_prop=0.3,
#'     nfold=5, alpha=NULL)
ml_model <- function(
        processed_se, char='none',
        transform=c('none', 'log10', 'square', 'cube'),
        ranking_method=c('p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet'),
        ml_method=c('Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost'),
        split_prop=0.3, nfold=10, alpha=NULL){
    ## check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if (is.null(char) | isFALSE(all(.check_char(processed_se, char, type='ml'))) ) {
        stop("Wrong char input, you can view the available char list by list_lipid_char function.")
    }
    if (is.null(ranking_method) | isFALSE(ranking_method %in% c('p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet')) ){
        stop("ranking_method must be one of 'p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet'.")
    } else if (ranking_method=='ElasticNet'){
        if (!is.numeric(alpha) | isFALSE(.check_numeric_range(alpha, 0, 1)) ) {
            stop("alpha must be a numeric value between 0 and 1.")
        }
    }
    if (is.null(ml_method) | isFALSE(ml_method %in% c('Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost')) ){
        stop("ml_method must be one of 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost'.")
    } else if (ml_method=='ElasticNet'){
        if (!is.numeric(alpha) | isFALSE(.check_numeric_range(alpha, 0, 1)) ) {
            stop("alpha must be a numeric value between 0 and 1.")
        }
    }
    if (!is.numeric(split_prop) | isFALSE(.check_numeric_range(split_prop, 0.1, 0.5)) ) {
        stop("split_prop must be a numeric value between 0.1 and 0.5.")
    }
    if (!is.numeric(nfold) | isFALSE(.check_numeric_range(nfold, 1, NULL)) ) {
        stop("nfold must be a numeric value >= 1.")
    }
    if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
        stop("transform must be one of 'none', 'log10', 'square', or 'cube'.")
    }
    ## process
    ml_data <- .ml_process(processed_se, char, transform=transform, type="transpose")

    n_ROC <- 300
    feature_num <- (ncol(ml_data)-1)
    sele_feature_num <- c(2, 3, 5, 10, 20, 50, 100)
    if(feature_num < 100){
        sele_feature_num <- c(
            sele_feature_num[sele_feature_num < feature_num], feature_num)
    }
    # Create cross-validation splits
    cv <- rsample::mc_cv(ml_data, prop=split_prop, times=nfold, strata=NULL)
    # Initialize result storage
    cv_results <- list(
        model=list(), cm=list(), ROC=list(), PR=list(), meanROC=list(), varimp=list())
    best_model <- list()
    cv_feature_save <- list()
    # Perform cross-validation
    for (fold in seq_len(nfold)) {
        message("Processing CV fold", fold)
        # Split data into train and test sets
        train_data <- rsample::assessment(cv$splits[[fold]])
        test_data <- rsample::analysis(cv$splits[[fold]])
        fold_results <- .ml_fold(
            train_data, test_data, ranking_method, ml_method, sele_feature_num,
            fold, n_ROC, alpha)
        # Store results
        ## cv result
        for (result_type in names(cv_results)) {
            cv_results[[result_type]] <- c(cv_results[[result_type]], fold_results[[result_type]])
        }
        ## model result
        if (length(best_model) == 0) {
            best_model <- fold_results$best_model
        } else {
            for (i in seq_along(fold_results$best_model$models)) {
                if (fold_results$best_model$scores[i] > best_model$scores[i]) {
                    best_model$models[[i]] <- fold_results$best_model$models[[i]]
                    best_model$features[[i]] <- fold_results$best_model$features[[i]]
                    best_model$scores[i] <- fold_results$best_model$scores[i]
                }
            }
        }
        cv_feature_save[[fold]] <- fold_results$feature_save
    }
    cv_feature_save <- lapply(cv_feature_save, .convert_feature_name, sele_feature_num)
    best_model$features <- .convert_feature_name(best_model$features, sele_feature_num)
    names(best_model$models) <- as.character(sele_feature_num)
    # Post-process results
    for (result_type in names(cv_results)) {
        cv_results[[result_type]] <- do.call(rbind, cv_results[[result_type]]) %>%
            dplyr::mutate(ranking_method=ranking_method) %>%
            dplyr::select(ranking_method, ml_method, cv_fold, feature_num, dplyr::everything())
    }
    # Special processing
    cv_results$meanROC <- cv_results$meanROC %>%
        dplyr::group_by(feature_num, threshold) %>%
        dplyr::summarise(
            sensitivity=mean(sensitivity, na.rm=TRUE), specificity=mean(specificity, na.rm=TRUE),
            precision=mean(precision, na.rm=TRUE), recall=mean(recall, na.rm=TRUE)) %>%
        dplyr::mutate(ranking_method=ranking_method, ml_method = ml_method, cv_fold='mean') %>%
        dplyr::select(ranking_method, ml_method, cv_fold, feature_num, dplyr::everything())
    cv_results$varimp %<>%
        dplyr::mutate(feature=.lipid_name_replace(cv_results$varimp$feature, type="revert"))
        # dplyr::mutate(feature=gsub('___', '\\|', feature)) %>%
        # dplyr::mutate(feature=gsub('_', ' ', gsub('__', ':', feature)))
    # final output
    ml_se <- SummarizedExperiment::SummarizedExperiment(
        assays=SummarizedExperiment::assay(processed_se),
        rowData=SummarizedExperiment::rowData(processed_se),
        colData=SummarizedExperiment::colData(processed_se),
        metadata=list(
            char=char, transform=transform, ranking_method=ranking_method,
            ml_method=ml_method, nfold=nfold, feature_option=sele_feature_num,
            model_result=cv_results$model, confusion_matrix=cv_results$cm, roc_result=cv_results$ROC,
            pr_result=cv_results$PR, mean_roc_result=cv_results$meanROC,
            feature_importance_result=cv_results$varimp, selected_features=cv_feature_save,
            best_model=best_model$models, best_model_feature=best_model$features))
    return(ml_se)
    }

## Helper functions
.ml_fold <- function(train_data, test_data, ranking_method, ml_method, sele_feature_num, fold, n_ROC, alpha) {
    fold_results <- list(
        model=list(), cm=list(), ROC=list(), PR=list(), meanROC=list(), varimp=list(),
        best_model=list(), feature_save=list() )
    for (feature_set in seq_along(sele_feature_num)) {
        n_features <- sele_feature_num[feature_set]
        # Select features
        selected_features <- .feature_selection(train_data, ranking_method, topN=n_features, alpha)
        fold_results$feature_save[[feature_set]] <- selected_features
        # Prepare data with selected features
        feature_cols <- c(TRUE, colnames(train_data)[-1] %in% selected_features)
        train_subset <- train_data[feature_cols]
        test_subset <- test_data[feature_cols]
        # Run ML model
        ml_result <- .ml_coreFun(train_subset, test_subset, ml_method, n_ROC, alpha)
        # Store results
        result_types <- c("model", "cm", "ROC", "PR", "meanROC", "varimp")
        for (i in seq_along(result_types)) {
            result <- ml_result[[i]] %>% dplyr::mutate(cv_fold=fold, feature_num=n_features)
            fold_results[[result_types[i]]][[feature_set]] <- result
        }
        # Update best model if necessary
        model_score <- ml_result$ROC_result$ROC_AUC[1]+ml_result$PR_result$PR_AUC[1]
        if (length(fold_results$best_model$scores) < feature_set ||
            model_score > fold_results$best_model$scores[feature_set]) {
            fold_results$best_model$models[[feature_set]] <- ml_result$model
            fold_results$best_model$features[[feature_set]] <- fold_results$feature_save[[feature_set]]
            fold_results$best_model$scores[feature_set] <- model_score
        }
    }
    return(fold_results)
}

.feature_selection <- function(data, ranking_method, topN, alpha){
    if (ranking_method %in% c('p_value', 'pvalue_FC')) {
        result <- data %>%
            dplyr::mutate(group=ifelse(group==0, 'group0', 'group1')) %>%
            tidyr::gather(key=feature, value=value, -group) %>%
            dplyr::group_by(group, feature) %>%
            dplyr::summarise(value=list(value)) %>%
            tidyr::spread(group, value) %>% dplyr::group_by(feature)
        if(ranking_method == 'p_value') {
            res <- result %>% dplyr::mutate(
                p_value=tryCatch(
                    stats::t.test(
                        unlist(group0), unlist(group1), var.equal=TRUE)$p.value,
                    error=function(e){NA})) %>% dplyr::select(feature, p_value)
            res %<>% dplyr::arrange(p_value) %>% .[seq_len(topN),1] %>% unlist()
        } else {
            res <- result %>% dplyr::mutate(
                pvalue_FC=-log10(tryCatch(
                    stats::t.test(unlist(group0), unlist(group1), var.equal=TRUE)$p.value,
                    error=function(e){NA}))*abs(
                        mean(unlist(group0), na.rm=TRUE)/mean(unlist(group1), na.rm=TRUE))) %>%
                dplyr::select(feature, pvalue_FC)
            res %<>% dplyr::arrange(dplyr::desc(pvalue_FC)) %>%
                .[seq_len(topN),1] %>% unlist()
        }
    }else if(ranking_method == 'ROC'){
        res <- data[-1] %>% tidyr::gather(key=feature, value=value) %>%
            dplyr::group_by(feature) %>% dplyr::summarise(value=list(value)) %>%
            dplyr::mutate(group=list(data[[1]])) %>% dplyr::group_by(feature) %>%
            dplyr::mutate(
                auc=tryCatch(
                    as.numeric(pROC::roc(unlist(group), unlist(value), quiet=TRUE)$auc),
                    error=function(e){NA})) %>%
            dplyr::select(feature, auc) %>%
            dplyr::arrange(dplyr::desc(auc)) %>% .[seq_len(topN),1,drop=TRUE]
    }else if(ranking_method == 'Random_forest'){
        model <- ranger::ranger(
            group ~ ., data=data, num.trees=500, importance='impurity', probability=TRUE)
        result <- model$variable.importance %>% sort(decreasing=TRUE) %>% names()
    }else if(ranking_method == 'SVM'){
        model <- e1071::svm(
            group~., kernel='linear', data=data, type='C-classification', probability=TRUE)
        result <- stats::coef(model)[-1] %>% abs() %>% sort(decreasing=TRUE) %>% names()
    }else if(ranking_method == 'Lasso'){
        result <- .fs_lasso_ridge_eNet(data, ranking_method, alpha=1, maxit=1000)
    }else if(ranking_method == 'Ridge'){
        result <- .fs_lasso_ridge_eNet(data, ranking_method, alpha=0, maxit=100000)
    }else if(ranking_method == 'ElasticNet'){
        result <- .fs_lasso_ridge_eNet(data, ranking_method, alpha=alpha, maxit=1000)
    }
    if (ranking_method %in% c('Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet')){
        res <- result[seq_len(topN)]
    }
    return(res)
}

.convert_feature_name <- function(feature_list, sele_feature_num){
    feature_list <- lapply(feature_list, function(x) {
        x <- .lipid_name_replace(x, type="revert")
        # x <- gsub('___', '|', x)
        # x <- gsub('__', ':', x)
        # x <- gsub('_', ' ', x)
        return(x)
    })
    names(feature_list) <- as.character(sele_feature_num)
    return(feature_list)
}

.fs_lasso_ridge_eNet <- function(data, ranking_method, alpha, maxit){
    #nlambda_value <- ifelse(type=="Ridge", 50, NULL)
    alpha_value <- ifelse(ranking_method=='ElasticNet', alpha, ifelse(ranking_method=="Ridge", 0, 1))
    if (ranking_method=="Ridge") {
        model <- glmnet::cv.glmnet(
            x =as.matrix(data[-1]), y=data[[1]], alpha=alpha_value, nlambda=50,
            family="binomial", nfolds=5, maxit=maxit)
    } else {
        model <- glmnet::cv.glmnet(
            x=as.matrix(data[-1]), y=data[[1]], alpha=alpha_value, family="binomial",
            nfolds=5, maxit=maxit)
    }
    coef <- stats::coef(model,s='lambda.min') %>%
        as.matrix() %>% as.data.frame() %>% .[-1, ,drop=FALSE]
    result <- coef[[1]]
    names(result) <- rownames(coef)
    result %<>% abs() %>% sort(decreasing=TRUE) %>% names()
    return(result)
}

.ml_coreFun <- function(train_data, test_data, ml_method, n_ROC, alpha){
    if(ml_method == 'Random_forest'){
        model <- ranger::ranger(
            group ~ ., data=train_data, num.trees=500, importance='impurity', probability=FALSE)
        pred_prob <- stats::predict(model, test_data[-1])$predictions
        varimp <- model$variable.importance
    }else if(ml_method == 'SVM'){
        model <- e1071::svm(
            group~., kernel='linear', data=train_data, type='C-classification', probability=TRUE)
        pred_prob <- stats::predict(model, test_data[-1], probability=TRUE) %>%
            attributes() %>% .$probabilities %>% .[,2]
        varimp <- stats::coef(model)[-1] %>% abs() %>% sort(decreasing=TRUE)
    }else if(ml_method == 'Lasso'){
        model <- glmnet::cv.glmnet(
            x=as.matrix(train_data[-1]), y=train_data[[1]], alpha=1,
            family="binomial", nfolds=5, maxit=1000)
        pred_prob <- stats::predict(
            model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]
        varimp <- stats::coef(model,s='lambda.min') %>% as.matrix() %>%
            as.data.frame()%>% .[-1, ,drop=FALSE]
    }else if(ml_method == 'Ridge'){
        model <- glmnet::cv.glmnet(
            x =as.matrix(train_data[-1]), y=train_data[[1]], alpha=0, family="binomial",
            nfolds=5, maxit=100000, nlambda=50)
        pred_prob <- stats::predict(
            model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]
        varimp <- stats::coef(model,s='lambda.min') %>%
            as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
    }else if(ml_method == 'ElasticNet'){
        model <- glmnet::cv.glmnet(
            x=as.matrix(train_data[-1]), y=train_data[[1]], alpha=alpha,
            family="binomial", nfolds=5, maxit=1000)
        pred_prob <- stats::predict(
            model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]
        varimp <- stats::coef(model,s='lambda.min') %>%
            as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
    }else if(ml_method=='xgboost'){
        cv_model <- xgboost::xgb.cv(
            data=as.matrix(train_data[-1]), nrounds=100,max_depth=3,
            label=train_data[[1]], nfold=5, verbose=FALSE, prediction=TRUE,
            objective="binary:logistic", early_stopping_rounds=10)
        cv_nround_min <- cv_model$best_iteration
        model <- xgboost::xgboost(
            data=as.matrix(train_data[-1]), label=train_data[[1]],
            nrounds=cv_nround_min, max_depth=3, objective="binary:logistic", verbose=0)
        pred_prob <- stats::predict(model, as.matrix(test_data[-1]))
        varimp <- xgboost::xgb.importance(colnames(test_data[-1]), model)
    }

    if (ml_method=='xgboost') {
        feature <- varimp$Feature
        importance <- varimp$Gain
    } else if (ml_method %in% c('Random_forest', 'SVM')) {
        feature <- names(varimp)
        importance <- varimp
    } else if (ml_method %in% c('Lasso', 'Ridge', 'ElasticNet')) {
        feature <- rownames(varimp)
        importance <- abs(varimp[[1]])
    }
    model_varImp <- .tab_model_varImp(ml_method, test_data, pred_prob, feature, importance)
    model_result <- model_varImp$model
    varImp_result <- model_varImp$varImp
    cm <- caret::confusionMatrix(
        as.factor(model_result$true_label),
        factor(model_result$pred_label, levels=c(0, 1)), positive='1')
    two_class_roc <- data.frame(
        truth=model_result$true_label, Class1=model_result$pred_prob) %>%
        dplyr::mutate(truth=ifelse(truth == 1, 'Class1', 'Class2')) %>%
        dplyr::mutate(truth=as.factor(truth))
    ROC_auc <- yardstick::roc_auc(two_class_roc, truth, Class1)$.estimate
    ROC_curve <- yardstick::roc_curve(two_class_roc, truth, Class1)
    PR_auc <- yardstick::pr_auc(two_class_roc, truth, Class1)$.estimate
    PR_curve <- yardstick::pr_curve(two_class_roc, truth, Class1)
    ROC_plot <- pROC::roc(model_result$true_label, model_result$pred_prob, quiet=TRUE)
    thershold <- c(seq(0, 1, length.out=n_ROC))
    meanROC <- pROC::coords(
        ROC_plot, thershold,
        ret=c("threshold", "specificity", "sensitivity", "precision", "recall"),
        transpose=FALSE)

    ## output
    cm_result <-rbind(data.frame(index=names(cm$overall)[1], value=cm$overall[1]),
                      data.frame(index=names(cm$byClass), value=cm$byClass)) %>%
        dplyr::mutate(value=value*100) %>% dplyr::mutate(ml_method=ml_method) %>%
        dplyr::select(ml_method, dplyr::everything())
    ROC_result <- data.frame(
        ml_method=ml_method, sensitivity=ROC_curve$sensitivity,
        specificity=ROC_curve$specificity, ROC_AUC=ROC_auc)
    PR_result <- data.frame(
        ml_method=ml_method, precision=PR_curve$precision,
        recall=PR_curve$recall, PR_AUC=PR_auc)
    meanROC_result <- meanROC %>% dplyr::mutate(ml_method=ml_method) %>%
        dplyr::select(ml_method, dplyr::everything())
    return(list(
        model_result=model_varImp$model, confusion_matrix=cm_result,
        ROC_result=ROC_result, PR_result=PR_result, meanROC_result=meanROC_result,
        varImp_result=model_varImp$varImp, model=model))
}

.tab_model_varImp <- function(ml_method, test_data, pred_prob, feature, importance){
    model <- data.frame(
        ml_method=ml_method, ID=rownames(test_data), true_label=test_data[[1]],
        pred_prob=pred_prob) %>% dplyr::mutate(pred_label=ifelse(pred_prob > 0.5,1,0))
    varImp <- data.frame(ml_method=ml_method, feature=feature, importance=importance)
    return(list(model=model, varImp=varImp))
}
