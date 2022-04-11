#' @title ROC_plot_all
#' @description The receiver operating characteristic (ROC) curve is one of the common methods to evaluate the diagnostic ability of a binary classifier. Mean AUC and 95% confidence interval for the ROC curve are calculated from all CV runs in each feature number.
#' @description This function provides
#' \enumerate{
#' \item the overall ROC Curve of CVs with different feature numbers.
#' \item the ROC Curve of average CVs by user-selected feature numbers.
#' }
#' @param data1 receiver operating characteristic after cross-validation. An output data frame of \code{\link{ML_final}} (list\bold{[[3]]}, cv_ROC_result).
#' @param data2 mean receiver operating characteristic of cross-validation. An output data frame of \code{\link{ML_final}} (list\bold{[[5]]}, cv_meanROC_result).
#' @param feature_n A numeric value specifying the number of features to be shown.
#' @return 1 tibble, 1 data frame, and 2 plot.
#' \enumerate{
#' \item data2[-c(8,9)]: tibble of ROC values
#' \item mean_AUC_plot: ROC curve plot
#' \item ROC_plot_data: ROC data frame of n features
#' \item ROC_plot: average ROC curve plot of n features
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
#' ROC_plot_all(ML_output[[3]], ML_output[[5]], feature_n=10)
ROC_plot_all <- function(data1,data2,feature_n){

  total_feature <- sort(unique(data1$feature_num)) %>% as.character()
  AUC_label <- c(length(total_feature))
  for(a in seq_len(length(total_feature))){
    cv_ROC_AUC <- data1[c(1,2,3,4,7)] %>% dplyr::filter(feature_num==total_feature[a]) %>%
      unique() %>% .$ROC_AUC
    cv_ROC_AUC <- stats::t.test(cv_ROC_AUC, mu=0.5)
    lower95 <- cv_ROC_AUC$conf.int[1] %>% round(3)
    upper95 <- cv_ROC_AUC$conf.int[2] %>% round(3)
    AUC_pvalue <- cv_ROC_AUC$p.value
    AUC_pvalue <- formatC(AUC_pvalue, digits = 2, format = "e")
    mean_AUC <- cv_ROC_AUC$estimate %>% round(3)
    AUC_label[a] <- stringr::str_c('AUC=', as.character(mean_AUC),
                          ' (',as.character(lower95), '-',
                          as.character(upper95),')')
  }
  label <- c(length(total_feature))

  for(a in seq_len(length(total_feature))){
    if(stringr::str_length(total_feature[a])==1){
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =5),' ', AUC_label[a])
    }else if(stringr::str_length(total_feature[a])==2){
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =4),' ', AUC_label[a])
    }else{
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =3),' ', AUC_label[a])
    }
  }

  mean_AUC_plot <- data2[c(3,4,6,7)] %>%
    dplyr::mutate(feature_num=as.character(feature_num))
  mean_AUC_plot$feature_num_label <- NA
  mean_AUC_plot$feature_num_label1 <- NA
  for(a in seq_len(length(total_feature))){
    mean_AUC_plot$feature_num_label[which(mean_AUC_plot$feature_num == total_feature[a])] <-label[a]
    mean_AUC_plot$feature_num_label1[which(mean_AUC_plot$feature_num == total_feature[a])] <-AUC_label[a]

  }
  #mean_AUC_plot$feature_num <- as.numeric(mean_AUC_plot$feature_num)
  mean_AUC_plot <- plotly::plot_ly(mean_AUC_plot,x=~(1-specificity),y=~sensitivity)%>%
    plotly::add_lines(name = ~reorder(feature_num_label,as.numeric(feature_num)),hoverinfo = "text",
              text = ~paste("feature_num :", feature_num,"<br>",feature_num_label1)) %>%
    plotly::add_lines(x = 0:1, y = 0:1, line = list(color = 'black', width = 2, dash = 'dash'),showlegend = FALSE)%>%
    plotly::layout(xaxis = list(title ="1-specificity"),
           legend = list(title=list(text="Feature number"),y=0.5,x=1.1,font = list(size = 9)),
           title =list(size =15,y=0.99,x=0.1,text="ROC curve"))

  ROC_plot_data <- data1 %>% dplyr::filter(feature_num==feature_n)
  cv_ROC_AUC <- data1[c(1,2,3,4,7)] %>% dplyr::filter(feature_num==feature_n) %>%
    unique() %>% .$ROC_AUC
  cv_ROC_AUC <- stats::t.test(cv_ROC_AUC, mu=0.5)
  lower95 <- cv_ROC_AUC$conf.int[1] %>% round(3)
  upper95 <- cv_ROC_AUC$conf.int[2] %>% round(3)
  AUC_pvalue <- cv_ROC_AUC$p.value
  AUC_pvalue <- formatC(AUC_pvalue, digits = 2, format = "e")
  mean_AUC <- cv_ROC_AUC$estimate %>% round(3)

  ROC_plot<- rbind(data1[3:6], data2[c(3,4,6,7)]) %>%
    dplyr::filter(feature_num==feature_n) %>%
    dplyr::mutate(cv=ifelse(cv_fold=='mean','1','0'))
  ROC_plot_1 <- ROC_plot %>% dplyr::filter(cv_fold!="mean")
  ROC_plot_2 <- ROC_plot %>% dplyr::filter(cv_fold=="mean")
  ROC_plot <- plotly::plot_ly(ROC_plot_1,x=~(1-specificity),y =~sensitivity,hoverinfo=NULL)%>%
    plotly::add_trace(data=ROC_plot_1,color=~cv_fold,colors ="gray",hoverinfo = "text",
              mode="lines",type='scatter',
              text = ~paste("cv_fold :", cv_fold))%>%
    plotly::add_trace(data=ROC_plot_2,name="mean",x=~(1-specificity),y =~sensitivity,line = list(color = "red"),
              hoverinfo = "text",mode="lines",type='scatter',
              text = ~paste("cv_fold :", cv_fold,"<br>AUC = ",
                            paste0(as.character(mean_AUC),'(',as.character(lower95),'-',as.character(upper95),')'))) %>%
    plotly::add_lines(x = 0:1, y = 0:1, line = list(color = 'black', width = 2, dash = 'dash'),showlegend = FALSE) %>%
    plotly::layout(xaxis = list(title ="1-specificity"),
           legend = list(title=list(text="cv_fold"),y=0.5,x=1.1,font = list(size = 9)),
           title =list(size =15,y=0.99,x=0.1,
                       text=stringr::str_c('ROC curve for ',as.character(feature_n),' feature model'))) %>%
    plotly::add_annotations(
      x=0.6,y=0.15,xref = "x",yref = "y",
      text = stringr::str_c('AUC=', as.character(mean_AUC),'(',as.character(lower95), '-',
                   as.character(upper95),')\n','pvalue=',as.character(AUC_pvalue)),
      xanchor = 'left',showarrow = FALSE)

  return(list(data2[-c(8,9)],mean_AUC_plot,ROC_plot_data,ROC_plot))
}

