#' @title Hclustering
#' @description Hierarchical clustering of lipid species or characteristics
#' derived from two groups and print a heatmap.
#' @param exp_data A data frame includes the expression of lipid features
#' in each sample. NAs are allowed. First column should be gene/lipid name and
#' first column name must be 'feature'.
#' @param DE_result_table A data frame comprises the significant results of
#' differential expression analysis, including fold change, p-value,
#' adjusted p-value. The output of \code{\link{DE_species_2}}.
#' @param group_info A data frame comprises the name of the sample, the label
#' of the sample, the group name of the sample, and the pair number represents
#' 'the pair' for the t-test/Wilcoxon test. NAs are allowed.
#' @param lipid_char_table A data frame. NAs are allowed. The name of first
#' column must be "feature".A data frame with lipid features, such as class,
#' total length. NAs are allowed. The name of first column must be "feature".
#' @param char_var A character string of the first lipid characteristic
#' selected by users from the column name of \bold{lipid_char_table},
#' such as total length.
#' @param distfun A character string of the distance measure indicating which
#' correlation coefficient (or covariance) is to be computed. Allowed methods
#' include \bold{"pearson"}, \bold{"kendall"}, and
#' \bold{"spearman"}.(default: "pearson")
#' @param hclustfun A character string of the agglomeration method to be used.
#' This should be (an unambiguous abbreviation of) one of \bold{"ward.D"},
#' \bold{"ward.D2"}, \bold{"single"}, \bold{"complete"},
#' \bold{"average"} (= UPGMA), \bold{"mcquitty"} (= WPGMA),
#' \bold{"median"} (= WPGMC), or
#' \bold{"centroid"} (= UPGMC). (default: "complete")
#' @param insert_ref_group A character string. The name of 'ctrl' after
#' name conversion.
#' @param ref_group A character string. The name of 'exp' after name conversion.
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @param type It should be 'all' or 'sig.' 'all' for output the results of all
#' lipid species or characteristics, and 'sig' for output the results of
#' significant lipid species or characteristics.
#' @return Return a list with a figures and a matrix.
#' \enumerate{
#' \item heatmap: a heatmap provides an overview of all/significant lipid
#' species or characteristics that illustrates the differences
#' between the control group and the experimental group.
#' \item data: the matrix of the heatmap.
#' }
#' @export
#' @examples
#' library(magrittr)
#' library(dplyr)
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' group_info <- DE_group_info
#' exp_transform <- data_process(exp_data, exclude_var_missing=TRUE,
#'                               missing_pct_limit=50, replace_zero=TRUE,
#'                               zero2what='min', xmin=0.5, replace_NA=TRUE,
#'                               NA2what='min', ymin=0.5, pct_transform=TRUE,
#'                               data_transform=TRUE, trans_type='log',
#'                               centering=FALSE,  scaling=FALSE)
#' exp_transform_non_log <- data_process(exp_data, exclude_var_missing=TRUE,
#'                                       missing_pct_limit=50,
#'                                       replace_zero=TRUE,
#'                                       zero2what='min', xmin=0.5,
#'                                       replace_NA=TRUE, NA2what='min',
#'                                       ymin=0.5, pct_transform=TRUE,
#'                                       data_transform=FALSE, trans_type='log',
#'                                       centering=FALSE, scaling=FALSE)
#' lipid_char_filter <- lipid_char_table %>%
#'    filter(feature %in% exp_transform$feature)
#' DE_species_table_sig <- DE_species_2(exp_transform_non_log,
#'                                      data_transform=TRUE,
#'                                      group_info=group_info, paired=FALSE,
#'                                      test='t.test',
#'                                      adjust_p_method='BH',
#'                                      sig_stat='p.adj',
#'                                      sig_pvalue=0.05,
#'                                      sig_FC=2)$DE_species_table_sig
#' char_var <- colnames(lipid_char_filter)[-1]
#' Hclustering(exp_transform, DE_result_table=DE_species_table_sig,
#'             group_info=group_info, lipid_char_table=lipid_char_filter,
#'             char_var=char_var[1], distfun='pearson',
#'             hclustfun='complete', plotly=TRUE, type='all')
Hclustering <- function(exp_data, DE_result_table, group_info,
                        lipid_char_table=NULL, char_var=NULL,
                        distfun='pearson', hclustfun='complete',
                        insert_ref_group=NULL, ref_group=NULL,
                        plotly=TRUE, type='all'){
  
  exp_data <- as.data.frame(exp_data)
  DE_result_table <- as.data.frame(DE_result_table)
  group_info <- as.data.frame(group_info)
  if(!is.null(lipid_char_table)){
    lipid_char_table <- as.data.frame(lipid_char_table)
  }
  if(ncol(exp_data) == 2){
    if(sum(class(exp_data[, -1]) %in% c("numeric", "integer")) != 1){
      stop("exp_data first column type must be 'character',
           others must be 'numeric'")
    }
  }else{
    if(sum(vapply(exp_data[, -1], class, character(1)) %in%
           c("numeric", "integer")) != ncol(exp_data[, -1])){
      stop("exp_data first column type must be 'character',
           others must be 'numeric'")
    }
  }
  if(tibble::is_tibble(exp_data)){
    if(nrow(exp_data) != nrow(unique(exp_data[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(exp_data) != length(unique(exp_data[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if(ncol(exp_data) < 3){
    stop("exp_data at least 2 samples.")
  }else if(ncol(exp_data) == 3){
    warning("exp_data only 2 samples will not show p-value,
            dotchart will color by log2FC")
  }
  if(nrow(exp_data) < 2){
    stop("exp_data number of lipids names (features) must be more than 2.")
  }
  if(sum(!is.na(exp_data[, -1])) == 0 | sum(!is.null(exp_data[, -1]))==0){
    stop("exp_data variables can not be all NULL/NA")
  }
  if(ncol(group_info) == 4){
    if(sum(vapply(group_info[, seq_len(3)], class,
                  character(1)) != "character") == 0){
      if("pair" %in% colnames(group_info)){
        if(which(colnames(group_info) == "pair") != 4){
          stop("group_info column must arrange in order of sample_name,
               label_name, group, pair(optional).")
        }
      }else{
        stop("group_info column must arrange in order of sample_name,
             label_name, group, pair(optional).")
      }
    }else{
      stop("group_info first 3 columns must be characters.")
    }
    if(sum(!is.na(group_info[, 4])) != 0 |
       sum(table(group_info[, 4]) != 2) != 0 &
       sum(is.na(group_info[, 4])) != 0){
      stop("group_info each pair must have a specific number, staring from
           1 to N. Cannot have NA, blank, or skip numbers.")
    }
    if(sum(group_info[, 1] %in% colnames(exp_data)) != nrow(group_info) |
       sum(group_info[, 1] %in% colnames(exp_data)) != ncol(exp_data[, -1])){
      stop("group_info 'sample_name' must same as the name of
           samples of exp_data")
    }
    if(length(unique(group_info[, 3])) == 2){
      if(sum(table(group_info[, 3]) >= 1) != 2){
        stop("group_info column 'group' only can have 2 groups,
             and >= 1 sample for each group.")
      }
    }else{
      stop("group_info column 'group' only can have 2 groups,
           and >= 1 sample for each group.")
    }
  }else if(ncol(group_info) == 3){
    if("pair" %in% colnames(group_info)){
      stop("group_info column must arrange in order of sample_name,
           label_name, group, pair(optional).")
    }
    if(sum(vapply(group_info,class,character(1)) != "character") != 0){
      stop("group_info first 3 columns must be characters.")
    }
    if(sum(group_info[, 1] %in% colnames(exp_data)) != nrow(group_info) |
       sum(group_info[, 1] %in% colnames(exp_data)) != ncol(exp_data[, -1])){
      stop("group_info 'sample_name' must same as the name of samples of
           exp_data")
    }
    if(length(unique(group_info[, 3])) == 2){
      if(sum(table(group_info[, 3]) >= 1) != 2){
        stop("group_info column 'group' only can have 2 groups,
             and >= 1 sample for each group.")
      }
    }else{
      stop("group_info column 'group' only can have 2 groups,
           and >= 1 sample for each group.")
    }
    if(!is.null(insert_ref_group)){
      if(!insert_ref_group %in% group_info[, 3]){
        stop("The insert_ref_group entered by users must be included
             in the group_info.")
      }
    }
  }
  if(!is.null(lipid_char_table)){
    if(nrow(lipid_char_table) == nrow(exp_data)){
      if(sum(lipid_char_table[, 1] %in% exp_data[, 1]) !=
         nrow(lipid_char_table)){
        stop("The lipids names (features) of lipid_char_table table must
             same as exp_data.")
      }
    }else{
      stop("The row number of lipid_char_table table must same as exp_data.")
    }
    if(!is(lipid_char_table[, 1], 'character')){
      stop("lipid_char_table first column must contain a list of
           lipids names (features).")
    }
    if(nrow(lipid_char_table) != length(unique(lipid_char_table[, 1]))){
      stop("lipid_char_table lipids names (features) must be unique.")
    }
    if("class" %in% colnames(lipid_char_table)){
      if(!is(lipid_char_table[, 'class'], 'character')){
        stop("lipid_char_table content of column 'class' must be characters")
      }
    }
    if("totallength" %in% colnames(lipid_char_table)){
      if(!class(lipid_char_table[, 'totallength']) %in%
         c("integer", "numeric")){
        stop("lipid_char_table content of column 'totallength' must be numeric")
      }
    }
    if("totaldb" %in% colnames(lipid_char_table)){
      if(!class(lipid_char_table[, 'totaldb']) %in% c("integer", "numeric")){
        stop("lipid_char_table content of column 'totaldb' must be numeric")
      }
    }
    if("totaloh" %in% colnames(lipid_char_table)){
      if(!class(lipid_char_table[, 'totaloh']) %in% c("integer", "numeric")){
        stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
      }
    }
    
    if(ncol(dplyr::select(lipid_char_table,
                          tidyselect::starts_with("FA_"))) != 0){
      FA_lipid_char_table <- lipid_char_table %>%
        dplyr::select(feature, tidyselect::starts_with("FA_"))
      FA_col <- grep("FA_", colnames(FA_lipid_char_table), value=TRUE)
      max_comma <- 0
      for(i in seq_len(length(FA_col))){
        col <- FA_col[i]
        comma_count <- max(stringr::str_count(FA_lipid_char_table[,col], ','),
                           na.rm=TRUE)
        if(comma_count > 0){
          FA_lipid_char_table <- FA_lipid_char_table %>%
            tidyr::separate(col,c(col,paste0(col, "_",
                                             seq_len(comma_count))),
                            ",", convert=TRUE)
        }
        if(comma_count > max_comma){max_comma <- comma_count}
      }
      FA_lipid_char_table <- FA_lipid_char_table %>%
        tidyr::gather(lipid.category, lipid.category.value,-feature)
      if(max_comma > 0){
        for (i in seq_len(max_comma)) {
          select_name <- paste0("_",i)
          FA_lipid_char_table <-FA_lipid_char_table[-intersect(
            grep(select_name, FA_lipid_char_table[, "lipid.category"]),
            which(is.na(FA_lipid_char_table$lipid.category.value))), ]
        }
      }
      if(is(FA_lipid_char_table$lipid.category.value, 'character') |
         sum(stats::na.omit(
           as.numeric(FA_lipid_char_table$lipid.category.value)) !=
           round(stats::na.omit(
             as.numeric(FA_lipid_char_table$lipid.category.value)))) != 0 |
         min(stats::na.omit(as.numeric(
           FA_lipid_char_table$lipid.category.value))) < 0){
        stop("In the 'FA_' related analyses, the values are positive integer or
             zero and separated by comma. i.e., 10,12,11")
      }
    }
  }
  
  colnames(exp_data)[1] <- 'feature'
  colnames(DE_result_table)[1] <- 'feature'
  
  rownames(exp_data) <- NULL
  exp.mat.all <- exp_data %>%
    dplyr::select(feature, group_info$sample_name) %>%
    tibble::column_to_rownames(var='feature') %>%
    as.matrix()
  colnames(exp.mat.all) <- group_info$label_name
  
  
  exp.mat.sig <- exp_data %>%
    dplyr::select(feature, group_info$sample_name) %>%
    dplyr::filter(feature %in% DE_result_table$feature) %>%
    tibble::column_to_rownames(var='feature') %>%
    as.matrix()
  colnames(exp.mat.sig) <- group_info$label_name
  
  if(!is.null(insert_ref_group) & !is.null(ref_group)){
    exp_raw_name <- ref_group[-which(insert_ref_group == ref_group)]
    group_info$group[which(group_info$group == 'ctrl')] <-  insert_ref_group
    group_info$group[which(group_info$group == 'exp')] <-  exp_raw_name
  }
  
  colGroup <- data.frame(Sample=group_info$group, stringsAsFactors=FALSE)
  col_color <- grDevices::rainbow(length(unique(colGroup$Sample)))
  col_color_label <- colGroup$Sample
  for(j in seq(unique(colGroup$Sample))){
    col_color_label[which(colGroup$Sample ==
                            unique(colGroup$Sample)[j])] <- col_color[j]
  }
  if(!is.null(lipid_char_table) & !is.null(char_var)){
    
    rowGroup.all <- exp.mat.all %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='feature') %>%
      dplyr::select(feature) %>%
      dplyr::left_join(lipid_char_table, by='feature') %>%
      dplyr::select(tidyselect::all_of(char_var))
    row_all_color <- grDevices::rainbow(length(unique(rowGroup.all$class)))
    row_all_color_label <- rowGroup.all$class
    for(j in seq(unique(rowGroup.all$class))){
      row_all_color_label[which(rowGroup.all$class ==
                                  unique(rowGroup.all$class)[j])] <-
        row_all_color[j]
    }
    
    rowGroup.sig <- exp.mat.sig %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='feature') %>%
      dplyr::select(feature) %>%
      dplyr::left_join(lipid_char_table, by='feature') %>%
      dplyr::select(tidyselect::all_of(char_var))
    row_sig_color <- grDevices::rainbow(length(unique(rowGroup.sig$class)))
    row_sig_color_label <- rowGroup.sig$class
    for(j in seq(unique(rowGroup.sig$class))){
      row_sig_color_label[which(rowGroup.sig$class ==
                                  unique(rowGroup.sig$class)[j])] <-
        row_sig_color[j]
    }
    if(type=='all'){
      row_color_label <- row_all_color_label
    }else if(type=='sig'){
      row_color_label <- row_sig_color_label
    }
  }else{
    row_color_label <- NULL
  }
  if(distfun %in% c("pearson", "kendall", "spearman")){
    dist_fun <- function(x) stats::as.dist(1-stats::cor(t(x), method=distfun))
    hclust_fun <- function(x) stats::hclust(x,  method=hclustfun)
  }else{
    dist_fun <- function(x) stats::dist(x, method=distfun)
    hclust_fun <- function(x) stats::hclust(x,  method=hclustfun)
  }
  heatmap_color_scale <- function(data){
    data <- round(data, 3)
    if(max(data) <= 0 & min(data) < 0){
      over_median <- min(data)/2
      if(max(data) < over_median){
        color <-  grDevices::colorRampPalette(c("#157AB5", "#92c5de"))(n=1000)
      }else{
        color_rank <- round(max(data)/(min(data))*1000)
        color_scale <- grDevices::colorRampPalette(c("#0571b0",
                                                     "#92c5de",
                                                     "white"))(n=1000)
        color <- color_scale[color_rank:1000]
      }
    }else if(min(data) >= 0 & max(data) > 0){
      over_median <- max(data)/2
      if(min(data) > over_median){
        color <-  grDevices::colorRampPalette(c("#f4a582", "#ca0020"))(n=1000)
      }else{
        color_rank <- round(min(data)/(max(data))*1000)
        color_scale <- grDevices::colorRampPalette(c("white",
                                                     "#f4a582",
                                                     "#ca0020"))(n=1000)
        color <- color_scale[color_rank:1000]
      }
    }
    return(color)
  }
  
  if(type=='all'){
    data <- exp.mat.all
  }else if(type=='sig'){
    data <- exp.mat.sig
  }
  
  if(nrow(data) >= 2 & ncol(data) >= 2 &
     sum(is.na(data)) == 0){
    
    if(max(nchar(rownames(data))) < 10){
      all_row_text_size <- 0.1
    }else if(max(nchar(rownames(data))) >= 10 &
             max(nchar(rownames(data))) < 20){
      all_row_text_size <- 0.2
    }else if(max(nchar(rownames(data))) >= 20 &
             max(nchar(rownames(data))) <30){
      all_row_text_size <- 0.3
    }else if(max(nchar(rownames(data))) >= 30 &
             max(nchar(rownames(data))) < 40){
      all_row_text_size <- 0.4
    }else {
      all_row_text_size <- 0.5
    }
    if(max(nchar(colnames(data))) < 10){
      all_col_text_size <- 0.1
    }else if(max(nchar(colnames(data))) >= 10 &
             max(nchar(colnames(data))) < 20){
      all_col_text_size <- 0.2
    }else if(max(nchar(colnames(data))) >= 20 &
             max(nchar(colnames(data))) < 30){
      all_col_text_size <- 0.3
    }else if(max(nchar(colnames(data))) >= 30 &
             max(nchar(colnames(data))) < 40){
      all_col_text_size <- 0.4
    }else {
      all_col_text_size <- 0.5
    }
    data <- sweep(data, 1, rowMeans(data, na.rm=TRUE))
    data <- sweep(data, 1, apply(data, 1,
                                 sd, na.rm=TRUE), "/")
    if(sum(is.na(data)) > 0){
      data <- data[-which(is.na(data[, 2])), ]
    }
    cb_grid <- iheatmapr::setup_colorbar_grid(y_length =0.6,
                                              x_start=1, y_start=0.4)
    if(distfun %in% c("pearson", "kendall", "spearman")){
      col_dend <- stats::hclust(
        stats::as.dist(1-stats::cor(data, method=distfun)),
        method=hclustfun)
      row_dend <- stats::hclust(
        stats::as.dist(1-stats::cor(t(data), method=distfun)),
        method=hclustfun)
    }else{
      col_dend <- stats::hclust(
        stats::dist(t(data), method=distfun),method=hclustfun)
      row_dend <- stats::hclust(
        stats::dist(data, method=distfun), method=hclustfun)
    }
    if(plotly == TRUE){
      if(min(data) >= 0 || max(data) <= 0){
        heatmap <- iheatmapr::iheatmap(
          data, colors=heatmap_color_scale(data),
          colorbar_grid=cb_grid, scale="rows")
      }else{
        heatmap <- iheatmapr::iheatmap(
          data, colorbar_grid=cb_grid, scale="rows")
      }
      heatmap <- heatmap %>%
        iheatmapr::add_col_annotation(annotation=colGroup,
                                      side="top",
                                      show_colorbar=FALSE)
      if(!is.null(lipid_char_table) & !is.null(char_var)){
        heatmap <- heatmap %>%
          iheatmapr::add_row_annotation(annotation=rowGroup.all,
                                        side="right",
                                        show_colorbar=FALSE)
      }
      heatmap <- heatmap %>%
        iheatmapr::add_col_dendro(col_dend, side="top", reorder =TRUE) %>%
        iheatmapr::add_row_dendro(row_dend, side="right", reorder =TRUE)
      if(ncol(data) <= 50){
        heatmap <- heatmap %>%
          iheatmapr::add_col_labels(side="bottom",
                                    size=all_col_text_size)
      }
      if(nrow(data)<=50){
        heatmap <- heatmap %>%
          iheatmapr::add_row_labels(side="left", size=all_row_text_size)
      }
    }else{
      if(min(data) >= 0 || max(data) <= 0){
        all_color_scale <- heatmap_color_scale(data)
      }else{
        all_color_scale <- grDevices::colorRampPalette(c("#92c5de",
                                                         "white",
                                                         "#f4a582"))(n=2500)
      }
      if(!is.null(lipid_char_table) & !is.null(char_var)){
        stats::heatmap(
          data, Rowv=TRUE, Colv=FALSE,
          dendrogram='both', trace="none",
          col=all_color_scale,
          distfun=dist_fun,
          hclustfun=hclust_fun,
          ColSideColors=col_color_label,
          RowSideColors=row_color_label,
          main=NULL,
          margins=c(8,8),
          lwid=c(1, 9),
          scale='none')
      }else{
        stats::heatmap(
          data, Rowv=TRUE, Colv=FALSE,
          dendrogram='both', trace="none",
          col=all_color_scale,
          distfun=dist_fun,
          hclustfun=hclust_fun,
          ColSideColors=col_color_label,
          main=NULL,
          margins=c(8,8),
          lwid=c(1, 9),
          scale='none')
      }
      
      heatmap <- grDevices::recordPlot()
      grDevices::dev.off()
    }
    reorder.data <- data[rev(row_dend$order), col_dend$order]
    
  }else{
    
    heatmap <- NULL
    reorder.data <- NULL
  }
  
  
  
  
  return(list(heatmap=heatmap, data=reorder.data))
} #function
