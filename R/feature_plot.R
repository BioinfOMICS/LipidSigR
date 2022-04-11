#' @title feature_plot
#' @description The function computes and ranks the contribution of each feature based on user-defined feature number, and visualizes the feature importance. In the ‘Algorithm-based’ part, when users choose a certain feature number, the selected frequency of top N features from all CV runs will be displayed. For a Linear SVM, Lasso, Ridge, or ElasticNet model, the importance of each feature depends on the absolute value of their coefficients in the algorithm, while Random Forest and XGBoost use built-in feature importance results.
#' @param data1 feature after cross-validation. An output list of \code{\link{ML_final}} (list[[6]], cv_feature_save).
#' @param data2 variable importance after cross-validation. An output data frame of \code{\link{ML_final}} (list[[7]],cv_varimp_result).
#' @param feature_n A numeric value specifying the number of elements to be shown.
#' @param nfold A numeric value. An optional parameter for setting the number of folds for cross-validation. (default=10)
#' @return Return a list of 1 data frame, 1 tibble, and 2 plots.
#' \enumerate{
#' \item feature_freq_data: a data frame of the selected frequency
#' \item feature_freq_plot: selected frequency plot
#' \item feature_imp_data: a tibble of feature importance
#' \item feature_imp_plot: feature importance plot.
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
#'                            trans_type='log', centering=FALSE, scaling=FALSE)
#' ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest',
#'                       ML_method='Random_forest', split_prop=0.3, nfold=10)
#' feature_plot(ML_output[[6]], ML_output[[7]], feature_n=10, nfold=10)
feature_plot <- function(data1, data2, feature_n, nfold=10){

  if(!class(feature_n)%in%c("numeric","integer")){
    stop('feature_n must be integer')
  }else{
    if(!feature_n %in% unique(data2$feature_num)){
      stop(paste0('feature_n must be one of ',paste(unique(data2$feature_num),collapse=", ")))
    }
  }
  if(!class(nfold)%in%c("numeric","integer")){
    stop('nfold must be integer')
  }
  feature_loc <- which(names(data1[[1]])==as.character(feature_n))
  cv_sele_feature <- data1 %>% purrr::map(function(x){x[[feature_loc]]}) %>%
    unlist %>% table()

  feature_freq_data <- data.frame(feature=names(cv_sele_feature),
                                  sele_freq=as.numeric(cv_sele_feature)/nfold) %>%
    dplyr::mutate(feature_num=feature_n) %>%
    dplyr::mutate(ranking_method=data2[1,1],
           ML_method=data2[1,2]) %>%
    dplyr::select(ranking_method,ML_method,
                  feature_num, dplyr::everything())


  if(nrow(feature_freq_data)>=10){
    plot_feature_num <- 10
  }else{
    plot_feature_num <- nrow(feature_freq_data)
  }
  #new
  feature_freq_data <- feature_freq_data %>% dplyr::arrange(dplyr::desc(sele_freq))
  #
  feature_freq_plot <- feature_freq_data[seq_len(plot_feature_num),] %>%
    plotly::plot_ly(x=~(sele_freq*100),y=~reorder(feature, sele_freq),type = "bar", orientation = "h",
            marker = list(color = ~-sele_freq,colorscale="Blues"),
            text = ~paste("Feature :", feature,"<br>Selected frequency :", sele_freq*100,"%"),hoverinfo = "text") %>%
    plotly::layout(xaxis =list(title="Selected frequency (%)",showgrid=TRUE,nticks=20,showline=TRUE,mirror='all'),
           yaxis =list(title=""))

  feature_imp_data <- data2 %>% dplyr::filter(feature_num==feature_n) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(importance=mean(importance, na.rm=TRUE)) %>%
    dplyr::select(-cv_fold) %>% unique() %>%
    dplyr:: mutate(ranking_method=data2[1,1],
           ML_method=data2[1,2]) %>%
    dplyr::select(ranking_method,ML_method,
                  feature_num, dplyr::everything())
  if(nrow(feature_imp_data)>=10){
    plot_feature_num <- 10
  }else{
    plot_feature_num <- nrow(feature_imp_data)
  }
  #new
  feature_imp_data <- feature_imp_data %>% dplyr::arrange(dplyr::desc(importance))
  #
  feature_imp_plot <- feature_imp_data[seq_len(plot_feature_num),] %>%
    plotly::plot_ly(x=~importance,y=~reorder(feature, importance),type = "bar", orientation = "h",
            marker = list(color = ~-importance,colorscale="Blues"),
            text = ~paste("Feature :", feature,"<br>Average importance :",round(importance,2)),hoverinfo = "text") %>%
    plotly::layout(xaxis =list(title='Average importance',showgrid=TRUE,nticks=20,showline=TRUE,mirror='all'),
           yaxis =list(title=""))

  feature_imp_data$importance <- round(feature_imp_data$importance, 3)

  return(list(feature_freq_data, feature_freq_plot,
              feature_imp_data,feature_imp_plot))
}
