#' @title model_for_net
#' @description Select feature importance from the model of
#' \code{\link{ML_final}}. Two methods can be applied, namely
#' 'Algorithm-based' and 'SHAP analysis', to rank and the feature importance.
#' By choosing a certain feature number, specifying the number of features to
#' be computed, and the average feature importance of the top 10 features from
#' all CV runs will be displayed.
#' @param data The output data frame(list[[2]]) of
#' \code{\link{ML_data_process}}.
#' @param ML_method A character string for the machine learning method to
#' be computed. Allowed methods include 'Random_forest', 'SVM', 'Lasso',
#' 'Ridge', 'ElasticNet', 'xgboost'.
#' @param varimp_method A character string indicating the method to input
#' feature importance. Allowed methods are
#' \bold{Algorithm-based} and \bold{SHAP}.
#' @param best_model The output list of \code{\link{ML_final}}.
#' @param best_model_feature The output list of \code{\link{ML_final}}.
#' @param feature_num A numeric value specifying the number of features to
#' be computed.
#' @param nsim A positive integer indicating the times of simulation.
#' @return Return 1 data frame.
#' \enumerate{
#' \item varImp_result: data frame of variable importance.
#' \emph{NOTE: This result can be further used for \code{\link{cor_network}}}
#' }
#' @export
#' @examples
#' library(dplyr)
#' data("ML_exp_data")
#' data("ML_lipid_char_table")
#' data("ML_condition_table")
#' condition_table <- ML_condition_table[85:144, ]
#' exp_data <- ML_exp_data[1:40, ] %>%
#'     select(feature, condition_table$sample_name)
#' lipid_char_table <- ML_lipid_char_table[1:40, ]
#' char_var <- colnames(lipid_char_table)[-1]
#' ML_data <- ML_data_process(exp_data, group_info = condition_table,
#'                            lipid_char_table, char_var[1],
#'                            exclude_var_missing=TRUE, missing_pct_limit=50,
#'                            replace_zero=TRUE, zero2what='min', xmin=0.5,
#'                            replace_NA=TRUE, NA2what='min', ymin=0.5,
#'                            pct_transform=TRUE, data_transform=TRUE,
#'                            trans_type='log', centering=FALSE, scaling=FALSE)
#' ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
#'                       ML_method='Random_forest', split_prop=0.3, nfold=3)
#' model_for_net(ML_data[[2]], ML_method='Random_forest',
#'               varimp_method='Algorithm-based', ML_output[[8]],
#'               ML_output[[9]], feature_num=10, nsim=5)
model_for_net <- function(data, ML_method, varimp_method, best_model,
                          best_model_feature,feature_num, nsim){

  data <- as.data.frame(data)
  if(!varimp_method %in% c('Algorithm-based', 'SHAP analysis')){
    stop('varimp_method must be one of the strings "Algorithm-based",
         or "SHAP analysis".')
  }
  if(!ML_method %in% c('Random_forest', 'SVM', 'Lasso', 'Ridge',
                       'ElasticNet','xgboost')){
    stop('ML_method must be one of the strings "Random_forest", "SVM",
         "Lasso", "Ridge", "ElasticNet", or "xgboost".')
  }
  if(!class(feature_num) %in% c("numeric", "integer")){
    stop('feature_num must be integer')
  }else{
    if(!feature_num %in% names(best_model)){
      stop(paste0('feature_num must be one of ',
                  paste(names(best_model), collapse=", ")))
    }
  }
  if(!class(nsim) %in% c("numeric", "integer")){
    stop('nsim must be integer')
  }

  if(varimp_method == 'Algorithm-based'){
    num <- which(as.numeric(names(best_model)) == feature_num)
    if(ML_method == 'Random_forest'){
      varimp <- best_model[[num]]$variable.importance
      varImp_result <- data.frame(ML_method='Random_forest',
                                  feature=names(varimp), importance=varimp)
    }else if(ML_method == 'SVM'){
      varimp <- stats::coef(best_model[[num]])[-1]
      varImp_result <- data.frame(ML_method='SVM', feature=names(varimp),
                                  importance=varimp)
    }else if(ML_method == 'Lasso'){
      varimp <- stats::coef(best_model[[num]],s = 'lambda.min') %>%
        as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method='Lasso',
                                  feature=rownames(varimp),
                                  importance=varimp[[1]])
    }else if(ML_method == 'Ridge'){
      varimp <- stats::coef(best_model[[num]],s = 'lambda.min') %>%
        as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method='Ridge',
                                  feature=rownames(varimp),
                                  importance=varimp[[1]])
    }else if(ML_method == 'ElasticNet'){
      varimp <- stats::coef(best_model[[num]],s = 'lambda.min') %>%
        as.matrix() %>% as.data.frame()%>% .[-1, ,drop=FALSE]
      varImp_result <- data.frame(ML_method='ElasticNet',
                                  feature=rownames(varimp),
                                  importance=varimp[[1]])
    }else if(ML_method == 'xgboost'){
      varimp <- xgboost::xgb.importance(best_model[[num]]$feature_names,
                                        best_model[[num]])
      varImp_result <- data.frame(ML_method='xgboost', feature=varimp$Feature,
                                  importance=varimp$Gain)
    }
  }else{
    SHAP <- function(data, best_model, best_model_feature,
                     ML_method, feature_n, nsim){

      std1 <- function(x){
        return ((x - min(x, na.rm=TRUE))/
                  (max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
      }


      which_model <- which(names(best_model) == as.character(feature_n))
      model <- best_model[[which_model]]

      best_model_feature <- colnames(data) %in%
        best_model_feature[[which_model]]
      best_model_feature[1] <- TRUE
      data <- data[best_model_feature]

      if(ML_method == 'xgboost'){

        shap_values <- SHAPforxgboost::shap.values(xgb_model=model,
                                                   X_train=as.matrix(data[-1]))
        shap_score <- shap_values$shap_score

        shap_long <- SHAPforxgboost::shap.prep(xgb_model=model,
                                               X_train=as.matrix(data[-1]),
                                               top_n=NULL)

      }else if(ML_method == 'Random_forest'){

        pred <- function(object, newdata) {
          stats::predict(object, data=newdata)$predictions
        } #RF

        shap_score <- fastshap::explain(
          model,
          X=data[-1],
          nsim=nsim,
          pred_wrapper=pred)
      }else if(ML_method == 'SVM'){

        pred <- function(object, newdata) {
          attributes(stats::predict(object, newdata=newdata,
                                    probability=TRUE))$probabilities[,2]
        } #SVM

        shap_score <- fastshap::explain(
          model,
          X=data[-1],
          nsim=nsim,
          pred_wrapper=pred)
      }else if(ML_method %in% c('Lasso', 'Ridge', 'ElasticNet')){

        pred <- function(object, newdata) {
          stats::predict(object, newx=as.matrix(newdata), type='response')[,1]
        } #Lasso, Ridge, ElasticNet

        shap_score <- fastshap::explain(
          model,
          X=data[-1],
          nsim=nsim,
          pred_wrapper=pred)
      }

      shap_long <- shap_score %>%
        tidyr::gather(key='variable', value='value') %>%
        dplyr::mutate(ID=rep(seq_len(nrow(shap_score)),
                             ncol(shap_score)), rfvalue=unlist(data[-1])) %>%
        dplyr::group_by(variable) %>%
        dplyr::mutate(stdfvalue=std1(rfvalue)) %>%
        dplyr::mutate(mean_value=mean(abs(value))) %>%
        dplyr::mutate(variable=as.factor(variable)) %>%
        dplyr::select(ID, variable, value, rfvalue, stdfvalue, mean_value) %>%
        data.table::as.data.table()

      shap_long <- as.data.frame(shap_long)
      colnames(shap_long)[c(3, 4, 5, 6)] <-
        c('shapley_value', 'raw_value',
          'normalized_value', 'mean_shapley_value')
      return(list(shap_score, shap_long))

    }
    varImp_result <- SHAP(data=data, best_model=best_model,
                          best_model_feature=best_model_feature,
                          ML_method=ML_method, feature_n=feature_num,
                          nsim=nsim)[[2]]

    varImp_dir <- varImp_result %>%
      dplyr::mutate(dir=shapley_value*raw_value) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(dir=sum(dir))
    varImp_result <- varImp_result[c(2, 6)] %>% unique() %>%
      dplyr::left_join(varImp_dir,by='variable') %>%
      dplyr::mutate(mean_shapley_value=ifelse(dir > 0,
                                              mean_shapley_value,
                                              -mean_shapley_value))

    varImp_result <- data.frame(ML_method=ML_method,
                                feature=varImp_result$variable,
                                importance=varImp_result$mean_shapley_value) %>%
      unique()

  }

  return(varImp_result)
}

