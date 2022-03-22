#' @title PR_plot_all
#' @description The Precision-Recall (PR) curve is one of the common methods to evaluate the diagnostic ability of a binary classifier. Mean AUC and 95% confidence interval for the PR curve are calculated from all CV runs in each feature number. PR curve is more sensitive to data with highly skewed datasets and offers a more informative view of an algorithm's performance.
#' @description This function provides
#' \enumerate{
#' \item the overall PR Curve of CVs with different feature numbers.
#' \item the PR Curve of average CVs by user-selected feature numbers.
#' }
#' @param data1 Precision-Recall after cross-validation. An output data frame of \code{\link{ML_final}} (list\bold{[[4]]}, cv_PR_result).
#' @param data2 mean receiver operating characteristic of cross-validation. An output data frame of \code{\link{ML_final}} (list\bold{[[5]]}, cv_meanROC_result).
#' @param feature_n A numeric value specifying the number of feature to be shown.
#' @return Return a list of 1 tibble, 1 data frame, and 2 plots.
#' \enumerate{
#' \item data2[-c(6,7)]: tibble of precision and recall values of each sample
#' \item mean_AUC_plot: PR curve plot
#' \item PR_plot_data: data frame of the AUC, recall, and precision of PR of n features.
#' \item PR_plot: average PR curve plot of n feature.
#' }
#' @export
#' @examples
#' \dontrun{
#' data("ML_exp_data")
#' data("ML_lipid_char_table")
#' data("ML_condition_table")
#' exp_data <- ML_exp_data
#' lipid_char_table <- ML_lipid_char_table
#' condition_table <- ML_condition_table
#' char_var <- colnames(lipid_char_table)[-1]
#' ML_data <- ML_data_process(exp_data, group_info = condition_table, lipid_char_table,
#'                            char_var[1], exclude_var_missing=T, missing_pct_limit=50,
#'                            replace_zero=T, zero2what='min', xmin=0.5, replace_NA=T,
#'                            NA2what='min', ymin=0.5, pct_transform=T, data_transform=T,
#'                            trans_type='log', centering=F, scaling=F)
#' ML_output <- ML_final(ML_data[[2]], ranking_method='Random_forest', ML_method='Random_forest',
#'                       split_prop=0.3, nfold=10)
#' PR_plot_all(ML_output[[4]], ML_output[[5]], feature_n=10)
#' }
PR_plot_all <- function(data1,data2,feature_n){

  total_feature <- sort(unique(data1$feature_num)) %>% as.character()
  AUC_label <- c(length(total_feature))
  for(a in 1:length(total_feature)){
    cv_PR_AUC <- data1[c(1,2,3,4,7)] %>% dplyr::filter(feature_num==total_feature[a]) %>%
      unique() %>% .$PR_AUC
    cv_PR_AUC <- stats::t.test(cv_PR_AUC, mu=0.5)
    lower95 <- cv_PR_AUC$conf.int[1] %>% round(3)
    upper95 <- cv_PR_AUC$conf.int[2] %>% round(3)
    AUC_pvalue <- cv_PR_AUC$p.value
    AUC_pvalue <- formatC(AUC_pvalue, digits = 2, format = "e")
    mean_AUC <- cv_PR_AUC$estimate %>% round(3)
    AUC_label[a] <- stringr::str_c('AUC=', as.character(mean_AUC),
                          ' (',as.character(lower95), '-',
                          as.character(upper95),')')
  }
  label <- c(length(total_feature))

  for(a in 1:length(total_feature)){
    if(stringr::str_length(total_feature[a])==1){
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =5),' ', AUC_label[a])
    }else if(stringr::str_length(total_feature[a])==2){
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =4),' ', AUC_label[a])
    }else{
      label[a] <- stringr::str_c(stringr::str_pad(total_feature[a],width =3),' ', AUC_label[a])
    }
  }

  mean_AUC_plot <- data2[c(3,4,8,9)] %>%
    dplyr::mutate(feature_num=as.character(feature_num))
  mean_AUC_plot$feature_num_label <- NA
  mean_AUC_plot$feature_num_label1 <- NA
  for(a in 1:length(total_feature)){
    mean_AUC_plot$feature_num_label[which(mean_AUC_plot$feature_num == total_feature[a])] <-label[a]
    mean_AUC_plot$feature_num_label1[which(mean_AUC_plot$feature_num == total_feature[a])] <-AUC_label[a]
  }
  mean_AUC_plot <- plotly::plot_ly(mean_AUC_plot,x=~recall,y=~precision)%>%
    plotly::add_lines(name = ~reorder(feature_num_label,as.numeric(feature_num)),hoverinfo = "text",
              text = ~paste("feature_num :", feature_num,"<br>",feature_num_label1)) %>%
    plotly::layout(legend = list(title=list(text="Feature number"),y=0.5,x=1.1,font = list(size = 9)),
           title =list(size =15,y=0.99,x=0.1,text="PR curve"))

  PR_plot_data <- data1 %>% dplyr::filter(feature_num==feature_n)
  cv_PR_AUC <- data1[c(1,2,3,4,7)] %>% dplyr::filter(feature_num==feature_n) %>%
    unique() %>% .$PR_AUC
  cv_PR_AUC <- stats::t.test(cv_PR_AUC, mu=0.5)
  lower95 <- cv_PR_AUC$conf.int[1] %>% round(3)
  upper95 <- cv_PR_AUC$conf.int[2] %>% round(3)
  AUC_pvalue <- cv_PR_AUC$p.value
  AUC_pvalue <- formatC(AUC_pvalue, digits = 2, format = "e")
  mean_AUC <- cv_PR_AUC$estimate %>% round(3)

  PR_plot <- rbind(data1[3:6], data2[c(3,4,8,9)]) %>%
    dplyr::filter(feature_num==feature_n) %>%
    dplyr::mutate(cv=ifelse(cv_fold=='mean','1','0'))
  PR_plot_1 <- PR_plot %>% dplyr::filter(cv_fold!="mean")
  PR_plot_2 <- PR_plot %>% dplyr::filter(cv_fold=="mean")
  PR_plot <- plotly::plot_ly(PR_plot_1,x=~recall,y =~precision,hoverinfo = NULL)%>%
    plotly::add_lines(color =~cv_fold,colors = "gray",hoverinfo = "text",
              text = ~paste("cv_fold :", cv_fold))%>%
    plotly::add_lines(data=PR_plot_2,name="mean",x=~recall,y =~precision,line = list(color = "red"),hoverinfo = "text",
              text = ~paste("cv_fold :", cv_fold,"<br>",
                            paste0(as.character(mean_AUC),'(',as.character(lower95),'-',as.character(upper95),')'))) %>%
    plotly::layout(legend = list(title=list(text="cv_fold"),y=0.5,x=1.1,font = list(size = 12)),
           title = list(size =15,y=0.99,x=0.1,
                        text=stringr::str_c('PR curve for ',as.character(feature_n) ,' feature model'))) %>%
    plotly::add_annotations(
      x=0.2,y=0.3,xref = "x",yref = "y",
      text = stringr::str_c('AUC=', as.character(mean_AUC),'(',as.character(lower95), '-',
                   as.character(upper95),')\n','pvalue=',as.character(AUC_pvalue)),
      xanchor = 'left',showarrow = F)

  return(list(data2[-c(6,7)],mean_AUC_plot,PR_plot_data,PR_plot))
}


