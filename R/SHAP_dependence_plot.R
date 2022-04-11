#' @title SHAP_dependence_plot
#' @description This function visualizes the SHAP values against the feature values for each variable.
#' @param shap_long The output data frame of \code{\link{SHAP}} (list[[2]]).
#' @param x A character string of the feature name selected by users from the "variable" column of \bold{shap_long}.
#' @param y A character string of the feature name selected by users from the "variable" column of \bold{shap_long}.
#' @param color_var A character string of the feature name selected by users from the "variable" column of \bold{shap_long}.
#' @return Return a list of 1 plot.
#' \enumerate{
#' \item shap_dependence_plot: SHAP dependence plot.
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
#' SHAP_output <- SHAP(ML_data[[2]], best_model=ML_output[[8]],
#'                     best_model_feature=ML_output[[9]],
#'                     ML_method='Random_forest', feature_n=10, nsim=5)
#' SHAP_dependence_plot(SHAP_output[[2]],x="C38.6.PC", y="C38.6.PC",
#'                      color_var="C38.6.PC")
SHAP_dependence_plot <- function(shap_long,x,y,color_var){

  if(!x %in%unique(shap_long$variable)){
    stop("x must be one of  shap_long 'variable'.")
  }
  if(!y %in%unique(shap_long$variable)){
    stop("y must be one of  shap_long 'variable'.")
  }
  if(!color_var %in%unique(shap_long$variable)){
    stop("color_var must be one of  shap_long 'variable'.")
  }
  shap_dependence_plot <- data.frame(a=(shap_long %>% dplyr::filter(variable==x) %>% .$raw_value),
                                     b=(shap_long %>% dplyr::filter(variable==y) %>% .$shapley_value),
                                     c=(shap_long %>% dplyr::filter(variable==color_var) %>% .$raw_value)) %>%
    plotly::plot_ly(x=~a,hoverinfo =NULL) %>%
    plotly::add_trace(y=~b,color = ~c,hoverinfo = "text",mode="markers",type='scatter',
              colors =~colorRampPalette(c("#E6E6FF" , "#0000FF"))(length(c)),
              text =~paste(x,":",a,"<br>Shapley_value :",round(b,3),
                           "<br>Feature value :",round(c,3))) %>%
    plotly::add_lines(y=~fitted(loess(b~a)),hoverinfo ="none",showlegend = FALSE) %>%
    plotly::colorbar(title = stringr::str_c(color_var,"\n(Feature value)")) %>%
    plotly::layout(title = 'Shap dependence plot',
           xaxis =list(title=x,nticks=20,showline=TRUE,mirror='all'),
           yaxis =list(title = stringr::str_c("Shapley value for ", y),zeroline = FALSE,zerolinewidth = 0))

  return(shap_dependence_plot)

}


