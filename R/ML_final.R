#' @title ML_final
#' @description This function constructs the machine learning model based on the output data processed by \code{\link{ML_data_process}}.
#' @param ML_data The output data frame(list[[2]]) of \code{\link{ML_data_process}}.
#' @param ranking_method A character string for the ranking method to be computed. Allowed methods include 'p_value', 'pvalue_FC', 'ROC', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet'.
#' @param ML_method A character string for the machine learning method to be computed. Allowed methods include 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet', 'xgboost'.
#' @param split_prop A numeric value between 0.1 to 0.5 indicates the proportion of data to be retained for modeling/analysis.
#' @param nfold A numeric value. The original dataset is randomly partitioned into N fold equal size subsamples.
#' @param alpha (optional) A numeric value between 0 and 1. Only need when selecting 'ElasticNet' as ML_method. 0 for Ridge and 1 for Lasso.
#' @return Return a list of 6 data frames and 3 lists.
#' \enumerate{
#' \item cv_model_result: data frame, an machine learning model after cross-validation.
#' \item cv_cm_result: data frame, confusion matrix after cross-validation.
#' \item cv_ROC_result: data frame, receiver operating characteristic after cross-validation.
#' \item cv_PR_result: data frame, Precision-Recall after cross-validation.
#' \item cv_meanROC_result: data frame, mean receiver operating characteristic of cross-validation.
#' \item cv_feature_save: list of length 10, feature after cross-validation.
#' \item cv_varimp_result: data frame, variable importance after cross-validation.
#' \item best_model: list of length 7.
#' \item best_model_feature: list of length 7.
#' }
#' @export
#' @examples
#' data("ML_exp_data")
#' data("ML_lipid_char_table")
#' data("ML_condition_table")
#' exp_data <- ML_exp_data
#' lipid_char_table <- ML_lipid_char_table
#' condition_table <- ML_condition_table
#' char_var <- colnames(lipid_char_table)[-1]
#' ML_data <- ML_data_process(exp_data, group_info = condition_table,
#'                            lipid_char_table, char_var[1],
#'                            exclude_var_missing=TRUE, missing_pct_limit=50,
#'                            replace_zero=TRUE, zero2what='min', xmin=0.5,
#'                            replace_NA=TRUE, NA2what='min', ymin=0.5,
#'                            pct_transform=TRUE, data_transform=TRUE,
#'                            trans_type='log', centering=FALSE,
#'                            scaling=FALSE)
#' ML_final(ML_data[[2]],ranking_method='Random_forest',
#'          ML_method='Random_forest', split_prop=0.3, nfold=10)
ML_final <- function(ML_data,ranking_method, ML_method,
                     split_prop, nfold, alpha=NULL){
  
  if(!ranking_method %in% c('p_value','pvalue_FC','ROC','Random_forest','SVM','Lasso','Ridge','ElasticNet')){
    stop('ranking_method must be one of the strings "p_value", "pvalue_FC", "ROC", "Random_forest", "SVM", "Lasso", "Ridge", or "ElasticNet".')
  }
  if(!ML_method %in% c('Random_forest','SVM','Lasso','Ridge','ElasticNet','xgboost')){
    stop('ML_method must be one of the strings "Random_forest", "SVM", "Lasso", "Ridge", "ElasticNet", or "xgboost".')
  }
  if(!class(nfold)%in%c("numeric","integer")){
    stop('nfold must be integer')
  }
  if(ML_method =='ElasticNet'){
    if(!is.null(alpha)){
      if(!class(alpha)%in%c("numeric","integer")){
        stop('alpha must be between 0 and 1')
      }else{
        if(alpha<0 | alpha>1){
          stop('alpha must be between 0 and 1')
        }
      }
    }else{
      stop('alpha must be between 0 and 1')
    }
  }

  n_ROC <- 300
  feature_selection <- function(data, ranking_method, topN){

    if(ranking_method=='p_value'){
      result <- data %>% dplyr::mutate(group=ifelse(group==0,'group0', 'group1')) %>%
        tidyr::gather(key = feature, value = value, -group) %>%
        dplyr::group_by(group, feature) %>%
        dplyr::summarise(value = list(value)) %>%
        tidyr::spread(group, value) %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(p_value = tryCatch(stats::t.test(unlist(group0), unlist(group1), var.equal = TRUE)$p.value,
                                  error=function(e){NA})) %>%
        dplyr::select(feature, p_value)

      result <- result %>% dplyr::arrange(p_value)  %>% .[seq_len(topN),1] %>% unlist()
      return(result)

    }else if(ranking_method=='pvalue_FC'){

      result <- data %>% dplyr::mutate(group=ifelse(group==0,'group0', 'group1')) %>%
        tidyr::gather(key = feature, value = value, -group) %>%
        dplyr::group_by(group, feature) %>%
        dplyr::summarise(value = list(value)) %>%
        tidyr::spread(group, value) %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(pvalue_FC = -log10(tryCatch(stats::t.test(unlist(group0), unlist(group1), var.equal = TRUE)$p.value,
                                           error=function(e){NA}))*abs(mean(unlist(group0),na.rm=TRUE)/mean(unlist(group1), na.rm=TRUE))) %>%
        dplyr::select(feature, pvalue_FC)


      result <- result %>% dplyr::arrange(dplyr::desc(pvalue_FC)) %>% .[seq_len(topN),1] %>% unlist()

      return(result)

    }else if(ranking_method=='ROC'){

      result <- data[-1] %>%
        tidyr::gather(key = feature, value = value) %>%
        dplyr::group_by(feature) %>%
        dplyr::summarise(value = list(value)) %>%
        dplyr::mutate(group=list(data[[1]])) %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(auc = tryCatch(as.numeric(pROC::roc(unlist(group), unlist(value),quiet=TRUE)$auc),
                              error=function(e){NA})) %>%
        dplyr::select(feature, auc) %>%
        dplyr::arrange(dplyr::desc(auc)) %>% .[seq_len(topN),1,drop=TRUE]

      return(result)

    }else if(ranking_method=='Random_forest'){
      model <- ranger::ranger(group ~ ., data = data, num.trees = 500, importance='impurity', probability = TRUE)
      result <- model$variable.importance %>%
        sort(decreasing = TRUE) %>% names()
      return(result[seq_len(topN)])
    }else if(ranking_method=='SVM'){
      model <- e1071::svm(group~., kernel='linear', data=data, type='C-classification', probability=TRUE)
      result <- stats::coef(model)[-1] %>% abs() %>% sort(decreasing = TRUE) %>% names()
      return(result[seq_len(topN)])
    }else if(ranking_method=='Lasso'){
      model <- glmnet::cv.glmnet(x =as.matrix(data[-1]),
                         y = data[[1]],
                         alpha = 1,  # lasso
                         family = "binomial",nfolds = 5,maxit=1000)

      coef <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]

      result <- coef[[1]]
      names(result) <- rownames(coef)
      result <- result %>% abs() %>% sort(decreasing = TRUE) %>% names()

      return(result[seq_len(topN)])
    }else if(ranking_method=='Ridge'){
      model <- glmnet::cv.glmnet(x =as.matrix(data[-1]),
                         y = data[[1]],
                         alpha = 0, nlambda=50,  # lasso
                         family = "binomial",nfolds = 5,maxit=100000)

      coef <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      result <- coef[[1]]
      names(result) <- rownames(coef)
      result <- result %>% abs() %>% sort(decreasing = TRUE) %>% names()
      return(result[seq_len(topN)])
    }else if(ranking_method=='ElasticNet'){
      model <- glmnet::cv.glmnet(x =as.matrix(data[-1]),
                         y = data[[1]],
                         alpha = alpha,  # lasso
                         family = "binomial",nfolds = 5,maxit=1000)

      coef <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      result <- coef[[1]]
      names(result) <- rownames(coef)
      result <- result %>% abs() %>% sort(decreasing = TRUE) %>% names()
      return(result[seq_len(topN)])
    }
  }

  ML_main <- function(train_data, test_data, ML_method){

    if(ML_method=='Random_forest'){
      model <- ranger::ranger(group ~ ., data = train_data, num.trees = 500, importance='impurity', probability = FALSE)

      pred_prob <- stats::predict(model, test_data[-1])$predictions

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- model$variable.importance
      varImp_result <- data.frame(ML_method=ML_method, feature=names(varimp),importance=varimp)
    }else if(ML_method=='SVM'){

      model <- e1071::svm(group~., kernel='linear', data=train_data, type='C-classification', probability=TRUE)

      pred_prob <- stats::predict(model, test_data[-1], probability = TRUE) %>% attributes() %>%
        .$probabilities %>% .[,2]

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- stats::coef(model)[-1] %>% abs() %>%
        sort(decreasing = TRUE)
      varImp_result <- data.frame(ML_method=ML_method, feature=names(varimp),importance=varimp)
    }else if(ML_method=='Lasso'){

      model <- glmnet::cv.glmnet(x =as.matrix(train_data[-1]),
                         y = train_data[[1]],
                         alpha = 1,  # lasso
                         family = "binomial",nfolds = 5,maxit=1000)


      pred_prob <- stats::predict(model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method=ML_method, feature=rownames(varimp),importance=abs(varimp[[1]]))
    }else if(ML_method=='Ridge'){
      model <- glmnet::cv.glmnet(x =as.matrix(train_data[-1]),
                         y = train_data[[1]],
                         alpha = 0,  # lasso
                         family = "binomial",nfolds = 5,maxit=100000,
                         nlambda=50)


      pred_prob <- stats::predict(model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method=ML_method, feature=rownames(varimp),importance=abs(varimp[[1]]))
    }else if(ML_method=='ElasticNet'){
      model <- glmnet::cv.glmnet(x =as.matrix(train_data[-1]),
                         y = train_data[[1]],
                         alpha = alpha,  # lasso
                         family = "binomial",nfolds = 5,maxit=1000)


      pred_prob <- stats::predict(model, s='lambda.min', newx=as.matrix(test_data[-1]), type="response")[,1]

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- stats::coef(model,s = 'lambda.min') %>% as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method=ML_method, feature=rownames(varimp),importance=abs(varimp[[1]]))
    }else if(ML_method=='xgboost'){

      cv_model <- xgboost::xgb.cv(data = as.matrix(train_data[-1]),
                         nrounds = 100,
                         max_depth = 3,
                         label = train_data[[1]],
                         nfold = 5,
                         verbose = FALSE,
                         prediction = TRUE,
                         objective = "binary:logistic",
                         early_stopping_rounds = 10)

      cv_nround_min <- cv_model$best_iteration

      model <- xgboost::xgboost(data = as.matrix(train_data[-1]),
                       label = train_data[[1]], nrounds = cv_nround_min,
                       max_depth = 3,
                       objective = "binary:logistic",
                       verbose=0)

      pred_prob <- stats::predict(model, as.matrix(test_data[-1]))

      model_result <- data.frame(ML_method=ML_method, ID=rownames(test_data), true_label=test_data[[1]],
                                 pred_prob=pred_prob) %>%
        dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

      varimp <- xgboost::xgb.importance(colnames(test_data[-1]),model)
      varImp_result <- data.frame(ML_method=ML_method, feature=varimp$Feature,importance=varimp$Gain)
    }

    cm <- caret::confusionMatrix(as.factor(model_result$true_label), factor(model_result$pred_label, levels = c(0,1)),positive = '1')
    cm_result <-rbind(data.frame(index=names(cm$overall)[1], value=cm$overall[1]),
                      data.frame(index=names(cm$byClass),
                                 value=cm$byClass)) %>%
      dplyr::mutate(value=value*100) %>%
      dplyr::mutate(ML_method=ML_method) %>%
      dplyr::select(ML_method, dplyr::everything())

    two_class_roc <- data.frame(truth=model_result$true_label,
                                Class1=model_result$pred_prob) %>%
      dplyr::mutate(truth=ifelse(truth==1, 'Class1', 'Class2')) %>%
      dplyr::mutate(truth=as.factor(truth))

    ROC_auc <- yardstick::roc_auc(two_class_roc, truth, Class1)[[1,3]]
    ROC_curve <- yardstick::roc_curve(two_class_roc, truth, Class1)
    PR_auc <- yardstick::pr_auc(two_class_roc, truth, Class1)[[1,3]]
    PR_curve <- yardstick::pr_curve(two_class_roc, truth, Class1)

    ROC_result <- data.frame(ML_method=ML_method,
                             sensitivity=ROC_curve$sensitivity,
                             specificity=ROC_curve$specificity,
                             ROC_AUC=ROC_auc)

    PR_result <- data.frame(ML_method=ML_method,
                            precision=PR_curve$precision,
                            recall=PR_curve$recall,
                            PR_AUC=PR_auc)


    ROC_plot <- pROC::roc(model_result$true_label, model_result$pred_prob, quiet = TRUE)
    thershold <- c(seq(0,1,length.out = n_ROC))

    meanROC <- pROC::coords(ROC_plot,thershold ,ret = c("threshold", "specificity", "sensitivity",
                                                  "precision", "recall"), transpose = FALSE)

    meanROC_result <- meanROC %>%
      dplyr::mutate(ML_method=ML_method) %>%
      dplyr::select(ML_method, dplyr::everything())

    return(list(model_result,cm_result, ROC_result, PR_result,meanROC_result,varImp_result,model))
  }


  ML <- function(data, ranking_method, ML_method, split_prop){
    feature_num <- (ncol(data)-1)
    sele_feature_num <- c(2,3,5,10,20,50,100)
    if(feature_num<100){
      sele_feature_num <- c(sele_feature_num[sele_feature_num<feature_num],
                            feature_num)
    }

    cv <- rsample::mc_cv(data, prop = split_prop, times = nfold, strata = NULL)

    num <- nrow(rsample::analysis(cv$splits[[1]]))
    num1 <- num+1
    num2 <- num+2

    cv_model_result <- data.frame(ML_method=rep(NA, 110000),
                                  ID=rep(NA, 110000),
                                  true_label=rep(NA, 110000),
                                  pred_prob=rep(NA, 110000),
                                  pred_label=rep(NA, 110000),
                                  cv_fold=rep(NA, 110000),
                                  feature_num=rep(NA, 110000))

    cv_cm_result <-  data.frame(ML_method=rep(NA, 3000),
                                index=rep(NA, 3000),
                                value=rep(NA, 3000),
                                cv_fold=rep(NA, 3000),
                                feature_num=rep(NA, 3000))

    cv_ROC_result<- data.frame(ML_method=rep(NA, 110000),
                               sensitivity=rep(NA, 110000),
                               specificity=rep(NA, 110000),
                               ROC_AUC=rep(NA, 110000),
                               cv_fold=rep(NA, 110000),
                               feature_num=rep(NA, 110000))

    cv_PR_result<- data.frame(ML_method=rep(NA, 110000),
                              precision=rep(NA, 110000),
                              recall=rep(NA, 110000),
                              PR_AUC=rep(NA, 110000),
                              cv_fold=rep(NA, 110000),
                              feature_num=rep(NA, 110000))
    cv_meanROC_result<- data.frame(ML_method=rep(NA, 110000),
                                   threshold=rep(NA, 110000),
                                   sensitivity=rep(NA, 110000),
                                   specificity=rep(NA, 110000),
                                   precision=rep(NA, 110000),
                                   recall=rep(NA, 110000),
                                   cv_fold=rep(NA, 110000),
                                   feature_num=rep(NA, 110000))

    cv_varimp_result<- data.frame(ML_method=rep(NA, 30000),
                                  feature=rep(NA, 30000),
                                  importance=rep(NA, 30000),
                                  cv_fold=rep(NA, 30000),
                                  feature_num=rep(NA, 30000))

    m <- 0
    best_model <- list(length(sele_feature_num))
    best_model_feature <- list(length(sele_feature_num))
    best_ROC_PR <- numeric(length(sele_feature_num))

    cv_feature_save <- list(nfold)
    for(a in seq_len(nfold)){
      methods::show(stringr::str_c('CV fold ', as.character(a), ' done'))
      train_data <- rsample::assessment(cv$splits[[a]])
      test_data <- rsample::analysis(cv$splits[[a]])

      feature_save <- list(length(sele_feature_num))
      for(b in seq_len(length(sele_feature_num))){
        sele_feature <- feature_selection(train_data, ranking_method, sele_feature_num[b])

        feature_save[[b]] <- sele_feature

        feature_loc <- colnames(train_data) %in% sele_feature
        feature_loc[1] <- TRUE

        ML_result <- ML_main(train_data[feature_loc],test_data[feature_loc],ML_method)
        model_result <- ML_result[[1]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])
        cm_result <- ML_result[[2]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])
        ROC_result <- ML_result[[3]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])

        PR_result <- ML_result[[4]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])

        meanROC_result <- ML_result[[5]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])

        varimp_result <- ML_result[[6]] %>%
          dplyr::mutate(cv_fold=a,feature_num=sele_feature_num[b])

        model_ROC_PR <- ML_result[[3]]$ROC_AUC[1]+ML_result[[4]]$PR_AUC[1]

        if(model_ROC_PR>best_ROC_PR[b]){
          best_model[[b]] <- ML_result[[7]]
          best_model_feature[[b]] <- sele_feature
          best_ROC_PR[b] <- model_ROC_PR
        }




        cv_model_result[(1+num*m):(num+num*m),] <- model_result
        cv_cm_result[(1+12*m):(12+12*m),] <- cm_result

        cv_ROC_result[(1+num2*m):(num2*m+nrow(ROC_result)),] <- ROC_result
        cv_PR_result[(1+num1*m):(num1*m+nrow(PR_result)),] <- PR_result
        cv_meanROC_result[(1+n_ROC*m):(n_ROC+n_ROC*m),] <- meanROC_result

        cv_varimp_result[(1+100*m):(100*m+nrow(varimp_result)),] <- varimp_result

        m <- m+1
      }
      names(feature_save) <- as.character(sele_feature_num)
      cv_feature_save[[a]] <- feature_save
    }

    cv_model_result <- cv_model_result %>% dplyr::mutate(ranking_method=ranking_method) %>%
      dplyr::select(ranking_method, ML_method, cv_fold, feature_num, dplyr::everything()) %>%
      dplyr::filter(!is.na(ML_method))
    cv_cm_result <- cv_cm_result %>% dplyr::mutate(ranking_method=ranking_method) %>%
      dplyr::select(ranking_method, ML_method, cv_fold, feature_num, dplyr::everything()) %>%
      dplyr::filter(!is.na(ML_method))

    cv_ROC_result <- cv_ROC_result %>% dplyr::mutate(ranking_method=ranking_method) %>%
      dplyr::select(ranking_method, ML_method, cv_fold, feature_num, dplyr::everything()) %>%
      dplyr::filter(!is.na(ML_method))

    cv_PR_result <- cv_PR_result %>% dplyr::mutate(ranking_method=ranking_method) %>%
      dplyr::select(ranking_method, ML_method, cv_fold, feature_num, dplyr::everything()) %>%
      dplyr::filter(!is.na(ML_method))

    cv_meanROC_result <- cv_meanROC_result %>%
      dplyr::filter(!is.na(threshold)) %>%
      dplyr::group_by(feature_num, threshold) %>%
      dplyr::summarise(sensitivity=mean(sensitivity,na.rm=TRUE),
                specificity=mean(specificity,na.rm=TRUE),
                precision=mean(precision,na.rm=TRUE),
                recall=mean(recall,na.rm=TRUE)) %>%
      dplyr::mutate(ranking_method=ranking_method,ML_method=ML_method, cv_fold='mean') %>%
      dplyr::select(ranking_method, ML_method, cv_fold,feature_num, dplyr::everything())

    cv_varimp_result <- cv_varimp_result %>% dplyr::mutate(ranking_method=ranking_method) %>%
      dplyr::select(ranking_method, ML_method, cv_fold, feature_num, dplyr::everything()) %>%
      dplyr::filter(!is.na(ML_method))

    names(best_model) <- as.character(sele_feature_num)
    names(best_model_feature) <- as.character(sele_feature_num)

    return(list(cv_model_result,cv_cm_result,
                cv_ROC_result,cv_PR_result,
                cv_meanROC_result,cv_feature_save,
                cv_varimp_result,best_model,
                best_model_feature))

  }


  result <- ML(ML_data, ranking_method, ML_method, split_prop = split_prop)

  return(result)
}






