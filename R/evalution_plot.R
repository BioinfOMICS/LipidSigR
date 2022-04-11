#' @title evalution_plot
#' @description This function provides useful indicators for users to evaluate model performance, such as the average value and 95% confidence interval of accuracy, sensitivity (recall), specificity, positive predictive value (precision), negative predictive value, F1 score, prevalence, detection rate, detection prevalence, balanced accuracy.
#' @param data confusion matrix after cross-validation. A output data frame of \code{\link{ML_final}} (list[[2]],cv_cm_result).
#' @param method A character string indicating which evaluation method to be used. Allowed methods include "Accuracy", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive Value", "Precision", "Recall", "F1", "Prevalence", "Detection Rate", "Detection Prevalence", "Balance Accuracy".
#' @return Return a list of 1 tibble and 1plot.
#' \enumerate{
#' \item evalution_data: a tibble of model evaluation information.
#' \item evalution_plot: model performance plot
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
#' ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest',
#'                       ML_method='Random_forest', split_prop=0.3, nfold=10)
#' evalution_plot(ML_output[[2]], method='Accuracy')
evalution_plot <- function(data, method){

  if(!method %in% c('Accuracy','Sensitivity','Specificity','Pos Pred Value','Neg Pred Value','Precision','Recall','F1','Prevalence','Detection Rate','Detection Prevalence', 'Balanced Accuracy')){
    stop('ranking_method must be one of the strings "Accuracy", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1", "Prevalence", "Detection Rate", "Detection Prevalence", or "Balanced Accuracy".')
  }
  evalution_data <- data %>% dplyr::filter(index==method) %>%
    dplyr::group_by(feature_num) %>%
    dplyr::mutate(mean_value=round(mean(value),2)) %>%
    dplyr::arrange(feature_num, cv_fold)
  evalution_plot <- data %>% dplyr::filter(index==method) %>%
    dplyr::group_by(feature_num) %>%
    dplyr::summarise(mean_value=round(mean(value),2),
              sd_value = round(sd(value),2)) %>%
    dplyr::mutate(highest=ifelse(mean_value==max(.$mean_value), 'red', 'black')) %>%
    ungroup() %>%
    plotly::plot_ly(x = ~feature_num, y =~mean_value, line = list(color='black'),
            marker = ~list(color=highest,size = 10),
            type='scatter',mode = 'lines',showlegend = FALSE,
            error_y = ~list(array = sd_value,color = '#000000'),
            text = ~paste("feature_num :", feature_num,"<br>mean :", mean_value,"<br>sd :", sd_value),hoverinfo = "text")%>%
    plotly::add_annotations(text = ~mean_value,xref = "x",yref = "y",showarrow = TRUE,
                    arrowhead = 4,arrowsize = .3,ax = -20,ay = 50) %>%
    plotly::layout(xaxis =list(title="Number of features",showgrid=TRUE,nticks=20,showline=TRUE,mirror='all'),
           yaxis =list(title="value (%)",range = c(0,100),nticks=15,mirror='all'),
           title =list(size =15,y=0.99,x=0.7,text= method))

  evalution_data$value <- round(evalution_data$value, 2)
  evalution_data$mean_value <- round(evalution_data$mean_value, 2)

  return(list(evalution_data,evalution_plot))
}


