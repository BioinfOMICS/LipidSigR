#' @title probability_plot
#' @description The function computes and visualizes the average predicted probabilities of each sample in testing data from all CV runs and allows users to explore those incorrect or uncertain labels. \\
#' @description Return two plots.
#' \enumerate{
#' \item the distribution of predicted probabilities in two reference labels
#' \item a confusion matrix composed of sample number and proportion
#' }
#' @param data an machine learning model after cross-validation. An output data frame of \code{\link{ML_final}} (list[[1]], cv_model_result).
#' @param feature_n A numeric value specifying the number of features to be computed.
#' @return Return a list of 2 tibbles and 2 plots.
#' \enumerate{
#' \item cm_data: a tibble of confusion matrix.
#' \item probability_plot: plot, the distribution of predicted probabilities in two reference labels.
#' \item cm_plot: plot. A confusion matrix composed of sample number and proportion.
#' \item data: a tibble of predicted probability and labels.
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
#' ML_output <- ML_final(ML_data[[2]],ranking_method='Random_forest', ML_method='Random_forest',
#'                       split_prop=0.3, nfold=10)
#' probability_plot(ML_output[[1]], feature_n=10)
#' }
probability_plot <- function(data, feature_n){

  if(!class(feature_n)%in%c("numeric","integer")){
    stop('feature_n must be integer')
  }else{
    if(!feature_n %in% unique(data$feature_num)){
      stop(paste0('feature_n must be one of ',paste(unique(data$feature_num),collapse=", ")))
    }
  }
  data <- data %>% dplyr::group_by(ID, feature_num) %>%
    dplyr::mutate(pred_prob=mean(pred_prob,na.rm=T)) %>%
    dplyr::select(-3,-8) %>% unique() %>%
    dplyr::mutate(pred_label=ifelse(pred_prob>0.5,1,0))

  probability_plot <-  data %>% dplyr::filter(feature_num==feature_n) %>%
    dplyr::mutate(true_label=as.factor(true_label)) %>%
    plotly::plot_ly(x=~true_label,y=~pred_prob, type = 'violin', box = list(visible = F),points=F,showlegend=F,
            color =~true_label,colors=c("#132B43", "#56B1F7")) %>%
    plotly::add_markers(x=~jitter(as.numeric(paste(true_label))),y=~pred_prob,
                text = ~paste("Actual group :", true_label,"<br>Predicted probability :", round(pred_prob,2)),hoverinfo = "text")%>%
    plotly::layout(xaxis =list(zeroline = F,title="Actual group"),
           yaxis =list(zeroline = F,title="Predicted probabilities"),
           title =list(size =15,y=0.99,x=0.1,text="Average sample probability in all CVs"))

  cm_data <-  data %>% dplyr::filter(feature_num==feature_n)

  cm <- caret::confusionMatrix(as.factor(cm_data$true_label), as.factor(cm_data$pred_label))

  cm_d <- as.data.frame(cm$table) %>%
    dplyr::group_by(Reference) %>%
    dplyr::mutate(pct=round(Freq/sum(Freq),2))

  cm_plot <- cm_d %>%
    dplyr::mutate(label=stringr::str_c(as.character(Freq),' (',as.character(pct),')')) %>%
    ggplot2::ggplot(ggplot2::aes(x = Reference , y =  Prediction, fill = pct))+
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradient(low = "white", high = "#08306b")+
    ggplot2::theme_bw()+
    ggplot2::labs(x='Actual group', y='Predicted group',title='Confusion matrix')+
    ggplot2::geom_text(ggplot2::aes(label = label,color = pct > 0.5))+
    ggplot2::scale_color_manual(guide = "none", values = c("black", "white"))+
    ggplot2::guides(fill="none")+
    ggplot2::coord_equal()
  cm_plot <- plotly::ggplotly(cm_plot)

  cm_data$pred_prob <- round(cm_data$pred_prob, 3)

  return(list(cm_data, probability_plot, cm_plot, data))
}



