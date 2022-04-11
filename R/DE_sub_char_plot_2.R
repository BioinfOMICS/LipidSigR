#' @title DE_sub_char_plot_2
#' @description This \code{\link{DE_sub_char_plot_2}} is to visualise the result of \code{\link{DE_sub_char_2}}.
#' \enumerate{
#' \item The subgroup analysis of the first user-selected characteristics.
#' \item Print a box plot, line plots (raw & sqrt scale), and bar plots (raw & sqrt scale).
#' }
#' @param DE_split_char_table_all A data frame. The output of \code{\link{DE_sub_char_2}}.
#' @param DE_split_char_index A data frame. The output of \code{\link{DE_sub_char_2}}.
#' @param group_info A data frame comprises the name of the sample, the label of the sample, the group name of the sample, and the pair number represents ‘the pair’ for the t-test/Wilcoxon test. NAs are allowed. NAs are allowed.
#' @param char_var A character string of the first lipid characteristic selected by users from the column name of \bold{lipid_char_table}, such as total length.
#' @param split_var A character string of the second lipid characteristic selected by users from \bold{lipid_char_table} for the subgroup analysis, such as class. \emph{NOTE: This parameter will be used to split data before entering main lipid characteristic 'char_var'.}
#' @param split_class A character string selected by users from the second user-selected lipid characteristic (\bold{split_var}) to visualise the specific plot for the selected category of that characteristic. For example, when user set split_var= 'class', and split_class='PC', this indicates the results will be plots of PC.
#' @param insert_ref_group A character string. The name of 'ctrl' after name conversion.
#' @param ref_group A character string. The name of 'exp' after name conversion.
#' @return Return a list with 5 figures.
#' \enumerate{
#' \item a bar plot of \bold{split_class}
#' \item a line plot of \bold{split_class}
#' \item a box plot of \bold{split_class}
#' \item a the bar plot of \bold{split_class} with the sqrt scale
#' \item a line plot of \bold{split_class} with the sqrt scale
#' }
#' @export
#' @examples
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' group_info <- DE_group_info
#' char_var <- colnames(lipid_char_table)[-1]
#' DE.sub.char.2 <- DE_sub_char_2(exp_data, data_transform=TRUE,
#'                                lipid_char_table = lipid_char_table,
#'                                split_var = char_var[1],
#'                                char_var = char_var[4],
#'                                group_info = group_info,paired = FALSE,
#'                                sig_pvalue=0.05, sig_FC=2,
#'                                exclude_var_missing=TRUE,
#'                                missing_pct_limit=50, replace_zero=TRUE,
#'                                zero2what='min', xmin=0.5,replace_NA=TRUE,
#'                                NA2what='min', ymin=0.5,
#'                                pct_transform=TRUE, trans_type='log',
#'                                centering=FALSE, scaling=FALSE)
#' char.class <- unique(DE.sub.char.2[[2]][1])
#' DE_sub_char_plot_2(DE.sub.char.2[[2]],
#'                    DE_split_char_index = DE.sub.char.2[[3]],
#'                    group_info = group_info, char_var = char_var[4],
#'                    split_var = char_var[1],
#'                    split_class = char.class[5,],
#'                    insert_ref_group=NULL, ref_group=NULL)
DE_sub_char_plot_2 <- function(DE_split_char_table_all, DE_split_char_index, group_info,
                               char_var='Category', split_var, split_class,
                               insert_ref_group=NULL, ref_group=NULL){

  #### Plot ####

  CTRL.RES <- DE_split_char_table_all %>%
    dplyr::select(seq_len(2), sig, mean_ctrl, sd_ctrl) %>%
    dplyr::mutate(Group = 'Ctrl')
  colnames(CTRL.RES) <- c('Split_category', 'Category', 'Significant', 'Mean', 'SD', 'Group')

  EXP.RES <- DE_split_char_table_all %>%
    dplyr::select(seq_len(2), sig, mean_exp, sd_exp) %>%
    dplyr::mutate(Group = 'Exp')
  colnames(EXP.RES) <- c('Split_category', 'Category', 'Significant', 'Mean', 'SD', 'Group')


  ## Fig.1 bar chart
  barTab <- data.table::rbindlist(l = list(CTRL.RES, EXP.RES), use.names = TRUE, fill = TRUE)
  if(sum(is.na(as.numeric(barTab$Category)))==0){
    barTab$Category <- as.factor(as.numeric(barTab$Category))
  }

  splitTab <- barTab %>% dplyr::filter(Split_category == split_class)
  splitTab <- splitTab %>% dplyr::group_by(Category) %>% dplyr::mutate(max_error_bar = max(Mean+SD)) %>% dplyr::ungroup()
  splitTab$post_hoc_pvalue = NA
  for(i in seq_len(nrow(splitTab))){
    post_hoc_pvalue_data <- DE_split_char_table_all[which(DE_split_char_table_all[,1]==split_class),]
    splitTab$post_hoc_pvalue[i] <- post_hoc_pvalue_data$post_hoc_pvalue[which(splitTab$Category[i]==post_hoc_pvalue_data[,2])]
  }
  if(!is.null(insert_ref_group) & !is.null(ref_group)){
    exp_raw_name <- ref_group[-which(insert_ref_group==ref_group)]
    splitTab$Group[which(splitTab$Group=='Ctrl')] <-  insert_ref_group
    splitTab$Group[which(splitTab$Group=='Exp')] <-  exp_raw_name
    splitTab$Group <- factor(splitTab$Group, levels = c(insert_ref_group, exp_raw_name))

  }
  splitTab_sig <- splitTab %>% dplyr::filter(Significant=='yes') %>% dplyr::mutate(pvalue_text = ifelse(post_hoc_pvalue<=0.001 ,"***",
                                                                                         ifelse(post_hoc_pvalue<=0.01 ,"**",
                                                                                                ifelse(post_hoc_pvalue<=0.05 ,"*",""))))

  barPlot <- ggplot2::ggplot(data=splitTab, ggplot2::aes(x=Category, y=Mean,fill=Group)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values=c('lightslateblue','sienna2')) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean, ymax=Mean+SD),color="gray39", width=.9,position=ggplot2::position_dodge()) +
    ggplot2::geom_text(data=splitTab_sig, ggplot2::aes(x=Category, y=max_error_bar+5 , label = pvalue_text),color="red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = char_var)
  barggplotly <- plotly::ggplotly(barPlot)
  for (i in seq_len(length(unique(barPlot$data$Group)))) {
    n <- length(unique(barPlot$data$Group))
    data <- barPlot$data[which(barPlot$data$Group==unique(barPlot$data$Group)[i]),]
    barggplotly$x$data[[i]]$text <- paste0("Category :",data$Category,
                                           "\nMean :",round(data$Mean,3),
                                           "\nSD :",round(data$SD,3),
                                           "\nGroup :",data$Group)
    barggplotly$x$data[[i+n]]$text <- paste0("Category :",data$Category,
                                             "\nMean :",round(data$Mean,3),
                                             "\nSD :",round(data$SD,3),
                                             "\nGroup :",data$Group)
  }
  for (i in seq_len(length(barggplotly$x$data))) {
    text = stringr::str_split(barggplotly$x$data[[i]]$hovertext,"<br />max_error_bar")
    hovertext =list()
    if(length(text)>0){
      for (j in seq_len(length(text))) {
        hovertext[[j]] <- paste(text[[j]][1],"Significant : YES")
      }
      barggplotly$x$data[[i]]$hovertext <-hovertext
    }
  }
  barPlot <- barggplotly

  ## Fig.1-1 bar chart sqrt scale

  barPlot_sqrt <- ggplot2::ggplot(data=splitTab, ggplot2::aes(x=Category, y=Mean,fill=Group)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values=c('lightslateblue','sienna2')) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean, ymax=Mean+SD),color="gray39", width=.9,position=ggplot2::position_dodge()) +
    ggplot2::geom_text(data=splitTab_sig, ggplot2::aes(x=Category, y=max_error_bar+5 , label = pvalue_text),color="red") +
    ggplot2::scale_y_sqrt() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = char_var)
  bar_ggplotly_sqrt <- plotly::ggplotly(barPlot_sqrt)
  for (i in seq_len(length(unique(barPlot_sqrt$data$Group)))) {
    n <- length(unique(barPlot_sqrt$data$Group))
    data <- barPlot_sqrt$data[which(barPlot_sqrt$data$Group==unique(barPlot_sqrt$data$Group)[i]),]
    bar_ggplotly_sqrt$x$data[[i]]$text <- paste0("Category :",data$Category,
                                                 "\nMean :",round(data$Mean,3),
                                                 "\nSD :",round(data$SD,3),
                                                 "\nGroup :",data$Group)
    bar_ggplotly_sqrt$x$data[[i+n]]$text <- paste0("Category :",data$Category,
                                                   "\nMean :",round(data$Mean,3),
                                                   "\nSD :",round(data$SD,3),
                                                   "\nGroup :",data$Group)
  }
  for (i in seq_len(length(bar_ggplotly_sqrt$x$data))){
    text = stringr::str_split(bar_ggplotly_sqrt$x$data[[i]]$hovertext,"<br />max_error_bar")
    hovertext =list()
    if(length(text)>0){
      for (j in seq_len(length(text))) {
        hovertext[[j]] <- paste(text[[j]][1],"Significant : YES")
      }
      bar_ggplotly_sqrt$x$data[[i]]$hovertext <-hovertext
    }
  }
  barPlot_sqrt <- bar_ggplotly_sqrt

  ## Fig.3 trend plot
  linePlot <- ggplot2::ggplot(data=splitTab, ggplot2::aes(x=Category, y=Mean,group=Group,color=Group)) +
    ggplot2::geom_line(stat="identity",position=ggplot2::position_dodge(0.05))+
    ggplot2::scale_color_manual(values=c('lightslateblue','sienna2')) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean-SD, ymax=Mean+SD),color="gray39",
                  position=ggplot2::position_dodge(0.05)) +
    ggplot2::geom_text(data=splitTab_sig,ggplot2::aes(x=Category, y=max_error_bar+5, label = pvalue_text),color="red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = char_var)
  lineggplotly <- plotly::ggplotly(linePlot)
  for(i in seq_len(length(lineggplotly$x$data))){
    if(!is.null(lineggplotly$x$data[i]$name)){
      lineggplotly$x$data[i]$name = gsub("\\(","",stringr::str_split(lineggplotly$x$data[i]$name,",")[[1]][1])
    }
  }
  for (i in seq_len(length(unique(linePlot$data$Group)))) {
    n <- length(unique(linePlot$data$Group))
    data <- linePlot$data[which(linePlot$data$Group==unique(linePlot$data$Group)[i]),]
    lineggplotly$x$data[[i]]$text <- paste0("Category :",data$Category,
                                            "\nMean :",round(data$Mean,3),
                                            "\nSD :",round(data$SD,3),
                                            "\nGroup :",data$Group)
    if(sum(grepl("\\*",lineggplotly$x$data[[i+n]]$text))==0){
      lineggplotly$x$data[[i+n]]$text <- paste0("Category :",data$Category,
                                                "\nMean :",round(data$Mean,3),
                                                "\nSD :",round(data$SD,3),
                                                "\nGroup :",data$Group)
    }
  }
  for (i in seq_len(length(lineggplotly$x$data))) {
    text = stringr::str_split(lineggplotly$x$data[[i]]$hovertext,"<br />max_error_bar")
    hovertext =list()
    if(length(text)>0){
      for (j in seq_len(length(text))) {
        hovertext[[j]] <- paste(text[[j]][1],"Significant : YES")
      }
      lineggplotly$x$data[[i]]$hovertext <-hovertext
    }
  }
  linePlot <- lineggplotly

  ## Fig.3_1 trend plot sqrt scale
  linePlot_sqrt <- ggplot2::ggplot(data=splitTab, ggplot2::aes(x=Category, y=Mean,group=Group,color=Group)) +
    ggplot2::geom_line(stat="identity",position=ggplot2::position_dodge(0.05))+
    ggplot2::scale_color_manual(values=c('lightslateblue','sienna2')) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Mean-SD, ymax=Mean+SD),color="gray39",
                  position=ggplot2::position_dodge(0.05)) +
    ggplot2::geom_text(data=splitTab_sig,ggplot2::aes(x=Category, y=max_error_bar+5, label = pvalue_text),color="red") +
    ggplot2::scale_y_sqrt() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = char_var)
  line_ggplotly_sqrt <- plotly::ggplotly(linePlot_sqrt)
  for(i in seq_len(length(line_ggplotly_sqrt$x$data))){
    if(!is.null(line_ggplotly_sqrt$x$data[i]$name)){
      line_ggplotly_sqrt$x$data[i]$name = gsub("\\(","",stringr::str_split(line_ggplotly_sqrt$x$data[i]$name,",")[[1]][1])
    }
  }
  for (i in seq_len(length(unique(linePlot_sqrt$data$Group)))) {
    n <- length(unique(linePlot_sqrt$data$Group))
    data <- linePlot_sqrt$data[which(linePlot_sqrt$data$Group==unique(linePlot_sqrt$data$Group)[i]),]
    line_ggplotly_sqrt$x$data[[i]]$text <- paste0("Category :",data$Category,
                                                  "\nMean :",round(data$Mean,3),
                                                  "\nSD :",round(data$SD,3),
                                                  "\nGroup :",data$Group)
    if(sum(grepl("\\*",line_ggplotly_sqrt$x$data[[i+n]]$text))==0){
      line_ggplotly_sqrt$x$data[[i+n]]$text <- paste0("Category :",data$Category,
                                                      "\nMean :",round(data$Mean,3),
                                                      "\nSD :",round(data$SD,3),
                                                      "\nGroup :",data$Group)
    }
  }
  for (i in seq_len(length(line_ggplotly_sqrt$x$data))) {
    text = stringr::str_split(line_ggplotly_sqrt$x$data[[i]]$hovertext,"<br />max_error_bar")
    hovertext =list()
    if(length(text)>0){
      for (j in seq_len(length(text))) {
        hovertext[[j]] <- paste(text[[j]][1],"Significant : YES")
      }
      line_ggplotly_sqrt$x$data[[i]]$hovertext <-hovertext
    }
  }
  linePlot_sqrt <- line_ggplotly_sqrt

  ## Fig.2 box plot
  colnames(DE_split_char_index)[1] <- 'Split_category'

  boxTab <- DE_split_char_index %>%
    dplyr::filter(Split_category == split_class) %>%
    dplyr::select(-Split_category) %>%
    tibble::column_to_rownames(var = char_var) %>%
    t() %>% as.data.frame() %>%
    merge(group_info, by.x = 0, by.y = 'sample_name')
  colnames(boxTab)[2] <- 'Category'
  if(!is.null(insert_ref_group) & !is.null(ref_group)){
    exp_raw_name <- ref_group[-which(insert_ref_group==ref_group)]
    boxTab$group[which(boxTab$group=='ctrl')] <-  insert_ref_group
    boxTab$group[which(boxTab$group=='exp')] <-  exp_raw_name
  }
  t.test.pvalue <- tryCatch({stats::t.test(Category ~ group, data = boxTab, var.equal = TRUE)["p.value"]},
                            warning = function(w) {NA},error = function(e){NA})
  if(!is.na(t.test.pvalue)){
    if(t.test.pvalue<=0.05){
      group_name <- c(unique(boxTab$group)[1],paste0(unique(boxTab$group)[1],"0"),unique(boxTab$group)[2])
      group_max <- max(boxTab$Category) %>%unique()
      group_min <- min(boxTab$Category) %>%unique()
      if(t.test.pvalue<=0.05 &t.test.pvalue>0.01){
        t.text <- c('', '*','')
      }else if (t.test.pvalue<=0.01 & t.test.pvalue>0.001){
        t.text <- c('', '**','')
      }else{
        t.text <- c('', '***','')
      }
      boxTab_1 <- boxTab %>%
        dplyr::filter(group==unique(boxTab$group)[1])
      boxTab_2 <- boxTab %>%
        dplyr::filter(group==unique(boxTab$group)[2])
      boxPlot <- plotly::plot_ly() %>%
        plotly::add_bars(x = group_name,
                 y = rep(group_max+0.15,3),
                 opacity=1,
                 showlegend = FALSE,
                 marker=list(line = list(color='rgba(0,0,0,0'),
                             color = 'rgba(0,0,0,0'),
                 textfont = list(color = 'red'),
                 text =  t.text ,
                 hoverinfo = 'none',
                 textposition = 'outside',
                 legendgroup = "1") %>%
        plotly::add_lines(x = c(rep(paste(unique(boxTab$group)[1]),2),rep(paste(unique(boxTab$group)[2]),2)),
                  y = c(group_max+0.1,group_max+0.15,group_max+0.15,group_max+0.1),
                  showlegend = FALSE,
                  line = list(color = 'black'),
                  legendgroup = "1",
                  hoverinfo = 'none') %>%
        plotly::add_boxplot(data = boxTab_1,x = ~group, y = ~Category,
                    color = I('lightslateblue'),
                    name = unique(boxTab$group)[1],
                    boxpoints = "all", jitter = 0.85,  pointpos = 0,
                    marker = list(size = 5, opacity = 0.8)) %>%
        plotly::add_boxplot(data = boxTab_2,x = ~group, y = ~Category,
                    color = I('sienna2'),
                    name = unique(boxTab$group)[2],
                    boxpoints = "all", jitter = 0.85,  pointpos = 0,
                    marker = list(size = 5, opacity = 0.8)) %>%
        plotly::layout(title = split_class,
               xaxis = list(title = 'Group',
                            tickmode = 'array',
                            tickvals =  c(unique(boxTab$group)[1], '', unique(boxTab$group)[2]),
                            ticktext = c(unique(boxTab$group)[1], '', unique(boxTab$group)[2]),
                            titlefont = list(size = 16),
                            tickfont = list(size = 14)),
               yaxis = list(title = paste0(char_var, ' index'),
                            titlefont = list(size = 16),
                            tickfont = list(size = 14),
                            range = c( group_min, group_max+0.5)
               ),
               legend = list(font = list(size = 14),
                             y = 0.5),
               margin = list(l=70, r=70, b=80, t = 60))
    }else{
      boxPlot <- plotly::plot_ly(data = boxTab,
                         x = ~group,
                         y = ~Category,
                         type = 'box',
                         color = ~group,
                         colors = c('lightslateblue','sienna2'),
                         boxpoints = 'all',
                         jitter = 0.85,
                         pointpos = 0,
                         marker = list(size = 5, opacity = 0.8)) %>%
        plotly::layout(title = split_class,
               xaxis = list(title = 'Group',
                            titlefont = list(size = 16),
                            tickfont = list(size = 14)),
               yaxis = list(title = paste0(char_var, ' index'),
                            titlefont = list(size = 16),
                            tickfont = list(size = 14)
               ),
               legend = list(font = list(size = 14),
                             y = 0.5),
               margin = list(l=70, r=70, b=80, t = 60))
    }
  }else{
    boxPlot <- plotly::plot_ly(data = boxTab,
                       x = ~group,
                       y = ~Category,
                       type = 'box',
                       color = ~group,
                       colors = c('lightslateblue','sienna2'),
                       boxpoints = 'all',
                       jitter = 0.85,
                       pointpos = 0,
                       marker = list(size = 5, opacity = 0.8)) %>%
      plotly::layout(title = split_class,
             xaxis = list(title = 'Group',
                          titlefont = list(size = 16),
                          tickfont = list(size = 14)),
             yaxis = list(title = paste0(char_var, ' index'),
                          titlefont = list(size = 16),
                          tickfont = list(size = 14)
             ),
             legend = list(font = list(size = 14),
                           y = 0.5),
             margin = list(l=70, r=70, b=80, t = 60))
  }

  return(list(barPlot, linePlot, boxPlot,barPlot_sqrt, linePlot_sqrt))

} #function
