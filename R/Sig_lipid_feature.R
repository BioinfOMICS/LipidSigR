#' @title Sig_lipid_feature
#' @description \code{\link{Sig_lipid_feature}} shows the significant lipid species based on different lipid characteristics and visualizes the difference between control and experimental groups by applying log2 Fold Change.
#' @param DE_species_table_sig A data frame comprises the significant results of differential expression analysis, including fold change, p-value, adjusted p-value. The output of \code{\link{DE_species_2}}.
#' @param lipid_char_table A data frame. NAs are allowed. The name of first column must be "feature".A data frame with lipid features, such as class, total length. NAs are allowed. The name of the first column must be "feature".
#' @param char_var A character string of the first lipid characteristic selected by users from the column name of \bold{lipid_char_table}, such as total length.
#' @param sig_FC Numeric. Significance of the fold-change. (default: 2)
#' @return Return a list with 3 figures.
#' \enumerate{
#' \item a bar chart shows the significant groups (values) with mean fold change over 2 in the selected characteristics by colors
#' \item a lollipop plot to compare multiple values simultaneously and it aligns the log2(fold change) of all significant groups(values) within the selected characteristics
#' \item a word cloud of the count of each group(value) of the selected characteristics
#' }
#' @export
#' @examples
#' \dontrun{
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' group_info <- DE_group_info
#' exp_transform <- data_process(exp_data, exclude_var_missing=T, missing_pct_limit=50,
#'                               replace_zero=T, zero2what='min', xmin=0.5, replace_NA=T,
#'                               NA2what='min', ymin=0.5, pct_transform=T, data_transform=T,
#'                               trans_type='log', centering=F,  scaling=F)
#' exp_transform_non_log <- data_process(exp_data, exclude_var_missing=F, missing_pct_limit=50,
#'                                       replace_zero=T, zero2what='min', xmin=0.5, replace_NA=T,
#'                                       NA2what='min', ymin=0.5, pct_transform=T,
#'                                       data_transform=F, trans_type='log',
#'                                       centering=F, scaling=F)
#' lipid_char_filter <- lipid_char_table %>% filter(feature %in% exp_transform_non_log$feature)
#' char_var <- colnames(lipid_char_filter)[-1]
#' DE_species_table_sig <- DE_species_2(exp_transform_non_log, data_transform = T,
#'                                      group_info = group_info, paired = F,
#'                                      test = 't.test', adjust_p_method = 'BH',
#'                                      sig_stat = 'p.adj', sig_pvalue = 0.05,
#'                                      sig_FC = 2)$DE_species_table_sig
#' Sig_lipid_feature(DE_species_table_sig, lipid_char_filter, char_var[1], sig_FC = 2)
#' }
Sig_lipid_feature <- function(DE_species_table_sig, lipid_char_table, char_var, sig_FC = 2){
  if(nrow(DE_species_table_sig)<1){
    stop("Less than 1 significant lipid features")
  }
  if(class(lipid_char_table[,1])!="character"){
    stop("lipid_char_table first column must contain a list of lipids names (features).")
  }
  if(nrow(lipid_char_table)!=length(unique(lipid_char_table[,1]))){
    stop("lipid_char_table lipids names (features) must be unique.")
  }
  if("class" %in%colnames(lipid_char_table) & class(lipid_char_table[,'class'])!="character"){
    stop("lipid_char_table content of column 'class' must be characters")
  }
  if("totallength" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totallength'])%in%c("integer","numeric")){
    stop("lipid_char_table content of column 'totallength' must be numeric")
  }
  if("totaldb" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totaldb'])%in%c("integer","numeric")){
    stop("lipid_char_table content of column 'totaldb' must be numeric")
  }
  if("totaloh" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totaloh'])%in%c("integer","numeric")){
    stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
  }
  if(!char_var %in% colnames(lipid_char_table)){
    stop("char_var must be included in the lipid_char_table.")
  }
  if(ncol(dplyr::select(lipid_char_table,tidyselect::starts_with("FA_")))==0){
    warning("(OPTIONAL) lipid_char_table does not contain column names starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>% dplyr::select(feature,tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_",colnames(FA_lipid_char_table),value = TRUE)
    max_comma <- 0
    for(i in 1:length(FA_col)){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[,col], ','),na.rm = T)
      if(comma_count>0){
        FA_lipid_char_table <- tidyr::separate(FA_lipid_char_table,col,c(col,paste0(col,"_",1:comma_count)),",", convert = TRUE)
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>% tidyr::gather(lipid.category, lipid.category.value,-feature)
    if(max_comma>0){
      for (i in 1:max_comma) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <-FA_lipid_char_table[-intersect(grep(select_name,FA_lipid_char_table[,"lipid.category"]),which(is.na(FA_lipid_char_table$lipid.category.value))),]
      }
    }
    if(class(FA_lipid_char_table$lipid.category.value)=="character"){
      stop("In the 'FA_' related analyses, the values are positive integer or zero and separated by comma. i.e., 10,12,11")
    }else if(sum(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value))!=round(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value))))!=0 | min(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)))<0){
      stop("In the 'FA_' related analyses, the values are positive integer or zero and separated by comma. i.e., 10,12,11")
    }
  }
  char <- lipid_char_table %>%
    dplyr::select(feature, tidyselect::all_of(char_var))
  colnames(char)[2] <- 'characteristic'

  plot.tab <- DE_species_table_sig %>%
    dplyr::left_join(char, by = 'feature') %>%
    dplyr::filter(!is.na(characteristic)) %>%
    dplyr::arrange(characteristic)
  plot.tab$characteristic <- factor(plot.tab$characteristic, sort(unique(plot.tab$characteristic)))


  #### bar plot ####
  sig_class.log2FC <- plot.tab %>%
    dplyr::group_by(characteristic)%>%
    dplyr::summarise(log2FC.mean=mean(log2FC),
              log2FC.sd=sd(log2FC),
              log2FC.direction=log2FC.mean/abs(log2FC.mean),
              significant='NO')%>%
    dplyr::mutate(log2FC.meansd=log2FC.mean+log2FC.sd*log2FC.direction,
           significant=replace(x = significant, list = abs(log2FC.mean) > log2(sig_FC), values = 'YES'))

  p.log2fc.sig <- ggplot2::ggplot(sig_class.log2FC,
                                  ggplot2::aes(x=stats::reorder(characteristic,-log2FC.mean), y=log2FC.mean, fill=significant)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(ggplot2::aes(min=log2FC.mean, max=log2FC.meansd)) +
    ggthemes::theme_hc() +
    ggplot2::scale_fill_manual(breaks=c('NO','YES'), values =c('#666666', '#FF6347'))+
    ggplot2::labs(fill=paste0('> ', sig_FC, ' FC'),
         y='log2(Fold Change)',
         x=char_var,
         title='Significant lipids') + # (without infinite fold change)
    ggplot2::theme(plot.title=ggplot2::element_text(size=16, hjust = 0.5),
          axis.title=ggplot2::element_text(size=14),
          axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1,size=12),#, angle = 30),
          axis.text.y=ggplot2::element_text(size=12))

  p.log2fc.sig$data$log2FC.mean <- round(p.log2fc.sig$data$log2FC.mean,3)
  p.log2fc.sig$data$log2FC.sd <- round(p.log2fc.sig$data$log2FC.sd,3)
  p.log2fc.sig$data$log2FC.meansd <- round(p.log2fc.sig$data$log2FC.meansd,3)
  in.log2fc.sig.ggplotly <- plotly::ggplotly(p.log2fc.sig)

  for(i in 1:length(in.log2fc.sig.ggplotly$x$data)){
    in.log2fc.sig.ggplotly$x$data[[i]]$text =gsub('(characteristic, -log2FC.mean)','',in.log2fc.sig.ggplotly$x$data[[i]]$text)
    in.log2fc.sig.ggplotly$x$data[[i]]$text =gsub('reorder',paste0(char_var,' :'),in.log2fc.sig.ggplotly$x$data[[i]]$text)
    in.log2fc.sig.ggplotly$x$data[[i]]$text =gsub('\\(\\): ','',in.log2fc.sig.ggplotly$x$data[[i]]$text)
    in.log2fc.sig.ggplotly$x$data[[i]]$text =gsub('log2FC.mean: ','mean(log2FC) :',in.log2fc.sig.ggplotly$x$data[[i]]$text)
    in.log2fc.sig.ggplotly$x$data[[i]]$text =gsub('significant: ','significant :',in.log2fc.sig.ggplotly$x$data[[i]]$text)
  }

  in.log2fc.sig <- in.log2fc.sig.ggplotly

  #### lolipop chart ####
  sig.dotchart <- plot.tab %>%
    dplyr::group_by(characteristic)%>%
    dplyr::mutate(log2FC.mean=mean(log2FC))
  sig.dotchart$characteristic <- forcats::fct_reorder(sig.dotchart$characteristic,sig.dotchart$log2FC.mean, min)


  p.sig.dotchart <- ggpubr::ggdotchart(sig.dotchart, combine = T,
                               x = "characteristic", y = "log2FC",
                               rotate = TRUE,
                               color = "white",
                               sorting = "none",
                               add = "segments",
                               dot.size = 3,
                               legend.title = char_var,
                               xlab = " ",
                               ylab = "log2(Fold change)",
                               legend = "right",
                               ggtheme = ggpubr::theme_pubr()) +
    ggplot2::geom_point(ggplot2::aes(text=paste("Characteristic :",characteristic,"<br>","log2FC : ",round(log2FC,2)),
                   color = characteristic,size=3)) +
    ggplot2::guides(size = "none")
  in.sig.dotchart <- plotly::ggplotly(p.sig.dotchart,tooltip = "text")


  #### word cloud ####
  wc.tab <- plot.tab %>%
    dplyr::group_by(characteristic) %>%
    dplyr::summarise(freqs = dplyr::n())

  wc <- hwordcloud::hwordcloud(text = wc.tab$characteristic,
                   size = wc.tab$freqs,
                   theme = "gridlight")


  return(list(barPlot = in.log2fc.sig,
              lolipop = in.sig.dotchart,
              word = wc))

} #function
