#' @title SHAP
#' @description This function uses the Shapley Additive exPlanations (SHAP) approach to rank and visualize the feature importance based on a user-defined feature number.
#' @param data The output data frame(list[[2]]) of \code{\link{ML_data_process}}.
#' @param best_model The ouput list of \code{\link{ML_data_process}} (list[[8]]).
#' @param best_model_feature The ouput list of \code{\link{ML_data_process}} (list[[9]]).
#' @param ML_method A character string for the machine learning method to be computed. Allowed methods include 'xgboost', 'Random_forest', 'SVM', 'Lasso', 'Ridge', 'ElasticNet'.
#' @param feature_n A numeric value specifying the number of features to be shown.
#' @param nsim A positive integer indicating the times of simulation.
#' @return Return a list of 1 table, 1 data frame, and 2 plots.
#' \enumerate{
#' \item shap_score: table of SHAP values.
#' \item shap_long: data frame, long-format SHAP data.
#' \item mean_shapley_plot: SHAP feature importance plot.
#' \item all_shapley_plot: SHAP value plot.
#' }
#' @export
#' @examples
#' \donttest{
#' data("ML_exp_data")
#' data("ML_lipid_char_table")
#' data("ML_condition_table")
#' exp_data <- ML_exp_data
#' lipid_char_table <- ML_lipid_char_table
#' condition_table <- ML_condition_table
#' char_var <- colnames(lipid_char_table)[-1]
#' ML_data <- ML_data_process(exp_data, group_info = condition_table, lipid_char_table,
#'                            char_var[1], exclude_var_missing=T, missing_pct_limit=50,
#'                            replace_zero=T, zero2what='min', xmin=0.5,replace_NA=T,
#'                            NA2what='min', ymin=0.5, pct_transform=T, data_transform=T,
#'                            trans_type='log', centering=F, scaling=F)
#' ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest', ML_method='Random_forest',
#'                       split_prop=0.3, nfold=10)
#' SHAP(ML_data[[2]], best_model=ML_output[[8]], best_model_feature=ML_output[[9]],
#'      ML_method='Random_forest', feature_n=10, nsim=5)
#' }
SHAP <- function(data, best_model, best_model_feature,  ML_method ,feature_n, nsim){

  if(!ML_method %in% c('Random_forest','SVM','Lasso','Ridge','ElasticNet','xgboost')){
    stop('ML_method must be one of the strings "Random_forest", "SVM", "Lasso", "Ridge", "ElasticNet", or "xgboost".')
  }
  if(!class(feature_n)%in%c("numeric","integer")){
    stop('feature_n must be integer')
  }else{
    if(!feature_n %in% names(best_model)){
      stop(paste0('feature_n must be one of ',paste(unique(data2$feature_num),collapse=", ")))
    }
  }
  if(!class(nsim)%in%c("numeric","integer")){
    stop('nsim must be integer')
  }

  std1 <- function(x){
    return ((x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }


  which_model <- which(names(best_model)==as.character(feature_n))
  model <- best_model[[which_model]]

  best_model_feature <- colnames(data)%in%best_model_feature[[which_model]]
  best_model_feature[1] <- T
  data <- data[best_model_feature]

  if(feature_n>10){topN <- 10}else{topN <- feature_n}

  if(ML_method=='xgboost'){

    shap_values <- SHAPforxgboost::shap.values(xgb_model = model, X_train = as.matrix(data[-1]))
    shap_score <- shap_values$shap_score

    shap_long <- SHAPforxgboost::shap.prep(xgb_model = model, X_train = as.matrix(data[-1]), top_n=NULL)

  }else if(ML_method=='Random_forest'){

    pred <- function(object, newdata) {
      stats::predict(object, data = newdata)$predictions
    } #RF

    shap_score <- fastshap::explain(
      model,
      X = data[-1],
      nsim = nsim,
      pred_wrapper = pred
    )
  }else if(ML_method=='SVM'){

    pred <- function(object, newdata) {
      attributes(stats::predict(object,newdata =newdata, probability=TRUE))$probabilities[,2]
    } #SVM

    shap_score <- fastshap::explain(
      model,
      X = data[-1],
      nsim = nsim,
      pred_wrapper = pred
    )
  }else if(ML_method %in% c('Lasso', 'Ridge', 'ElasticNet')){

    pred <- function(object, newdata) {
      stats::predict(object, newx = as.matrix(newdata), type = 'response')[,1]
    } #Lasso, Ridge, ElasticNet

    shap_score <- fastshap::explain(
      model,
      X = data[-1],
      nsim = nsim,
      pred_wrapper = pred
    )
  }

  shap_long <- shap_score %>%
    tidyr::gather(key='variable', value='value') %>%
    dplyr::mutate(ID=rep(1:nrow(shap_score),ncol(shap_score)), rfvalue=unlist(data[-1])) %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(stdfvalue=std1(rfvalue)) %>%
    dplyr::mutate(mean_value=mean(abs(value))) %>%
    dplyr::mutate(variable=as.factor(variable)) %>%
    dplyr::select(ID,variable, value, rfvalue,stdfvalue, mean_value) %>%
    data.table::as.data.table()

  mean_shapley_plot <- shap_long[,c(2,6)] %>% as.data.frame %>% unique() %>%
    dplyr::arrange(dplyr::desc(mean_value)) %>% .[1:topN,]%>%
    plotly::plot_ly(x=~mean_value,y=~reorder(variable,mean_value),type = "bar",orientation = "h",
            hoverinfo = "text",marker = list(color = ~-mean_value,colorscale="Blues"),
            text =~paste("Feature :", variable,"<br>Mean shapley value :",round(mean_value,2)))%>%
    plotly::layout(title = "SHAP feature importance",
           xaxis =list(title="mean(|Shapley value|)",nticks=15,showline=TRUE,mirror='all'),
           yaxis =list(title=" "))

  all_shapley_plot <- SHAPforxgboost::shap.plot.summary.wrap2(shap_score = shap_score,X = as.matrix(data[-1]), top_n = topN)
  all_shapley_plot <- plotly::ggplotly(all_shapley_plot)
  shap_long <- as.data.frame(shap_long)
  colnames(shap_long)[c(3,4,5,6)] <- c('shapley_value', 'raw_value', 'normalized_value', 'mean_shapley_value')


  return(list(shap_score, shap_long, mean_shapley_plot,all_shapley_plot))

}

