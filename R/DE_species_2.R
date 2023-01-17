#' @title DE_species_2
#' @description Compute differentially expressed analysis of
#' two groups (independent) based on the input Group Information to find
#' significant lipid species. To find differentially expressed lipid species,
#' a statistical method (t-test or Wilcoxon test) and the cut-offs for
#' significant lipid species need to be chosen, and the p-value will then be
#' adjusted by Benjamini-Hochberg or other methods.
#' @param exp_data A data frame includes the expression of lipid features in
#' each sample. NAs are allowed. First column should be gene/lipid name and
#' first column name must be 'feature'.
#' @param data_transform Logical. If data_transform = TRUE,
#' transform exp_data by log10.
#' @param group_info A data frame comprises the name of the sample,
#' the label of the sample, the group name of the sample, and the pair number
#' represents 'the pair' for the t-test/Wilcoxon test. NAs are allowed.
#' @param paired Logical. If paired = TRUE,
#' data are paired samples. (default: FALSE)
#' @param test A character string indicating which method to be used for
#' comparing means. Allowed method include \bold{"t.test"} and
#' \bold{"wilcox.test"}. (default: "t.test")
#' @param adjust_p_method Correction method, a character string. One of
#' \bold{"holm"}, \bold{"hochberg"}, \bold{"hommel"}, \bold{"bonferroni"},
#' \bold{"BH"}, \bold{"BY"}, \bold{"fdr"}, \bold{"none"}, can be
#' abbreviated. (default: "BH")
#' @param sig_stat A character string indicating which p-value is to be used
#' for the statistically significant. One of \bold{"p.adj"} or
#' \bold{"p"}. (default: "p.adj")
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @param sig_FC Numeric. Significance of the fold-change. (default: 2)
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @return Return a list with 2 data frames and 1 figures.
#' \enumerate{
#' \item DE_species_table_all: a data frame comprises the results of
#' differential expression analysis, including fold change, p-value,
#' adjusted p-value.
#' \item DE_species_table_sig: a data frame comprises the significant
#' results of differential expression analysis, including fold change,
#' p-value, adjusted p-value.
#' \item DE_species_dotchart_sig: a lollipop chart reveals the lipid species
#' that pass chosen cut-offs.
#' \item DE_species_maplot: an MA plot of lipid species
#' \item DE_species_volcano: a volcano plot of lipid species
#' }
#' @export
#' @examples
#' data("DE_exp_data")
#' data("DE_lipid_char_table")
#' data("DE_group_info")
#' exp_data <- DE_exp_data
#' lipid_char_table <- DE_lipid_char_table
#' group_info <- DE_group_info
#' exp_transform_non_log <- data_process(exp_data, exclude_var_missing=TRUE,
#'                                       missing_pct_limit=50,
#'                                       replace_zero=TRUE, zero2what='min',
#'                                       xmin=0.5, replace_NA=TRUE,
#'                                       NA2what='min', ymin=0.5,
#'                                       pct_transform=TRUE,
#'                                       data_transform=FALSE,
#'                                       trans_type='log',
#'                                       centering=FALSE, scaling=FALSE)
#' DE_species_2(exp_transform_non_log, data_transform=TRUE,
#'              group_info=group_info, paired=FALSE,
#'              test='t.test', adjust_p_method='BH', sig_stat='p.adj',
#'              sig_pvalue=0.05, sig_FC=2, plotly=TRUE)
DE_species_2 <- function(exp_data, data_transform=TRUE, group_info,
                         paired=FALSE, test='t.test', adjust_p_method='BH',
                         sig_stat='p.adj', sig_pvalue=0.05, sig_FC=2,
                         plotly=TRUE){

  exp_data <- as.data.frame(exp_data)
  if(!is(exp_data[, 1], 'character')){
    stop("exp_data first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data) == 2){
    if(!is(exp_data[, 1], 'character') |
       sum(class(exp_data[, -1]) %in% c("numeric", "integer")) != 1){
      stop("exp_data first column type must be 'character',
           others must be 'numeric'")
    }
  }else{
    if(!is(exp_data[, 1], 'character') |
       sum(vapply(exp_data[, -1], class, character(1)) %in%
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
  if(nrow(exp_data) < 5){
    stop("exp_data number of lipids names (features) must be more than 5.")
  }
  if(sum(!is.na(exp_data[, -1])) == 0 | sum(!is.null(exp_data[, -1])) == 0){
    stop("exp_data variables can not be all NULL/NA")
  }
  if(ncol(group_info) == 4){
    if(sum(vapply(group_info[, seq_len(3)],
                  class, character(1)) != "character") == 0){
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
      stop("group_info each pair must have a specific number, staring from 1
           to N. Cannot have NA, blank, or skip numbers.")
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
    if(sum(vapply(group_info, class, character(1)) != "character") != 0){
      stop("group_info first 3 columns must be characters.")
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
  }

  colnames(exp_data)[1] <- 'feature'

  if(paired == TRUE){
    exp_data_ga <- exp_data %>%
      tidyr::gather(sample_name, value, -1) %>%
      dplyr::left_join(group_info, by='sample_name') %>%
      dplyr::arrange(group, pair, feature)
  }else if(paired == FALSE){
    exp_data_ga <- exp_data %>%
      tidyr::gather(sample_name, value, -1) %>%
      dplyr::left_join(group_info, by='sample_name')
  }


  if(data_transform == TRUE & test == 't.test'){

    exp_data_trans_ga <- exp_data_ga
    exp_data_trans_ga[,3] <- log10(exp_data_trans_ga[, 3])

    exp_data_tab <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
             mean_exp=mean(unlist(exp), na.rm=TRUE),
             method='t-test',
             FC=mean_exp/mean_ctrl,
             log2FC=log2(FC))

    exp_data_p <- exp_data_trans_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(p_value=tryCatch(
        stats::t.test(unlist(ctrl), unlist(exp),
                      paired=paired, var.equal=TRUE)$p.value,
        error=function(e){NA}),
             m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value,
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)

    exp_data_stat <- exp_data_tab[-2:-3] %>%
      dplyr::right_join(exp_data_p[-2:-3], by='feature') %>%
      dplyr::mutate(sig_p=ifelse(p_value < sig_pvalue &
                                   abs(log2FC) > log2(sig_FC), 'yes', 'no'),
             sig_p_adj=ifelse(p_adj < sig_pvalue &
                                abs(log2FC) > log2(sig_FC), 'yes', 'no'))


  }else if(data_transform == FALSE & test == 't.test'){

    exp_data_p <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
             mean_exp=mean(unlist(exp), na.rm=TRUE),
             method='t-test',
             FC=mean_exp/mean_ctrl,
             log2FC=log2(FC),
             p_value=tryCatch(
               stats::t.test(unlist(ctrl), unlist(exp),
                             paired=paired, var.equal=TRUE)$p.value,
               error=function(e){NA}),
             m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value,
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)

    exp_data_stat <- exp_data_p[-2:-3] %>%
      dplyr::mutate(sig_p=ifelse(p_value < sig_pvalue &
                                   abs(log2FC) > log2(sig_FC), 'yes', 'no'),
             sig_p_adj=ifelse(p_adj < sig_pvalue &
                                abs(log2FC) > log2(sig_FC), 'yes', 'no'))

  }else if(data_transform == TRUE & test == 'wilcox.test'){

    exp_data_trans_ga <- exp_data_ga
    exp_data_trans_ga[, 3] <- log10(exp_data_trans_ga[, 3])

    exp_data_tab <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
             mean_exp=mean(unlist(exp), na.rm=TRUE),
             method='wilcoxon test',
             FC=mean_exp/mean_ctrl,
             log2FC=log2(FC))

    exp_data_p <- exp_data_trans_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(p_value=tryCatch(
        stats::wilcox.test(unlist(ctrl), unlist(exp), paired=paired)$p.value,
        error=function(e){NA}),
        m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value,
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)

    exp_data_stat <- exp_data_tab[-2:-3] %>%
      dplyr::right_join(exp_data_p[-2:-3], by='feature') %>%
      dplyr::mutate(sig_p=ifelse(p_value < sig_pvalue &
                                   abs(log2FC) > log2(sig_FC), 'yes', 'no'),
             sig_p_adj=ifelse(p_adj < sig_pvalue &
                                abs(log2FC) > log2(sig_FC), 'yes', 'no'))


  }else if(data_transform == FALSE & test == 'wilcox.test'){

    exp_data_p <- exp_data_ga %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(value=list(value)) %>%
      tidyr::spread(group, value) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(mean_ctrl=mean(unlist(ctrl), na.rm=TRUE),
             mean_exp=mean(unlist(exp), na.rm=TRUE),
             method='wilcoxon test',
             FC=mean_exp/mean_ctrl,
             log2FC=log2(FC),
             p_value=tryCatch(
               stats::wilcox.test(unlist(ctrl), unlist(exp),
                                  paired=paired)$p.value,
               error=function(e){NA}),
             m_log10_p_value=-log10(p_value))
    exp_data_p$p_adj <- stats::p.adjust(exp_data_p$p_value,
                                        method=adjust_p_method)
    exp_data_p$m_log10_p_adj <- -log10(exp_data_p$p_adj)

    exp_data_stat <- exp_data_p[-2:-3] %>%
      dplyr::mutate(sig_p=ifelse(p_value < sig_pvalue &
                                   abs(log2(FC)) > log2(sig_FC), 'yes', 'no'),
             sig_p_adj=ifelse(p_adj < sig_pvalue &
                                abs(log2(FC)) > log2(sig_FC), 'yes', 'no'))

  }

  exp_data_stat$mean_ctrl <- round(exp_data_stat$mean_ctrl, 5)
  exp_data_stat$mean_exp <- round(exp_data_stat$mean_exp, 5)
  exp_data_stat$FC <- round(exp_data_stat$FC, 3)
  exp_data_stat$log2FC <- round(exp_data_stat$log2FC, 3)
  exp_data_stat$m_log10_p_value <- round(exp_data_stat$m_log10_p_value, 3)
  exp_data_stat$m_log10_p_adj <- round(exp_data_stat$m_log10_p_adj, 3)

  #### Significant lipid ####
  if(sum(!is.na(exp_data_stat[,'p_value'])) == 0){
    sig.diff.exp <- exp_data_stat %>%
      dplyr::filter(abs(log2FC) >= log2(sig_FC)) %>%
      dplyr::distinct(feature, .keep_all=TRUE)
    warning('Only two samples could not be tested, significant lipids were
            instead defined by log2 fold change.')
  }else{
    if(sig_stat == 'p'){
      sig.diff.exp <- exp_data_stat %>%
        dplyr::filter(sig_p == 'yes') %>%
        dplyr::distinct(feature, .keep_all=TRUE)
    }else if(sig_stat == 'p.adj'){
      sig.diff.exp <- exp_data_stat %>%
        dplyr::filter(sig_p_adj == 'yes') %>%
        dplyr::distinct(feature, .keep_all=TRUE)
    }
  }

  if(nrow(sig.diff.exp) > 0){

    #### lolipop chart ####
    if(sum(is.finite(sig.diff.exp$log2FC)) == 0){ ##all sig. lipid is infinite

      #sig.diff.infinite <- sig.diff.exp
      sig.diff.exp$log2FC <- ifelse(sig.diff.exp$log2FC > 0, 5, -5)

      p.dotchart <- ggpubr::ggdotchart(sig.diff.exp, x="feature", y="log2FC",
                               rotate=TRUE,
                               color="white",
                               sorting="descending",
                               add="segments",
                               dot.size=2.5,
                               legend.title="log2(FC)",
                               xlab=" ",
                               ylab="log2(Fold Change)",
                               legend="right",
                               ggtheme=ggpubr::theme_pubr()) #+
      p.dotchart <- if(sum(!is.na(sig.diff.exp[,'p_value'])) == 0){
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature,
                                    "<br>", "log2(FC) : ", round(log2FC, 2),
                                    "<br>", "-log10(p-value) :",
                                    m_log10_p_value), color=log2FC, size=2.5)) +
          ggplot2::guides(size='none') +
          ggplot2::labs(colour="log2(FC)") +
          ggplot2::scale_colour_gradient2(low="steelblue",
                                          mid="white", high="red", midpoint=0) +
          ggplot2::scale_y_continuous(breaks=c(-5, -2, -1, 0, 1, 2, 5),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-6, 6))
      }else if(sig_stat == "p"){
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature, "<br>",
                                    "log2(FC) : ", round(log2FC, 2), "<br>",
                                    "-log10(p-value) :", m_log10_p_value),
                         color=m_log10_p_value, size=2.5)) +
          ggplot2::guides(size="none") +
          ggplot2::labs(colour="-log10(p-value)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0) +
          ggplot2::scale_y_continuous(breaks=c(-5, -2, -1, 0, 1, 2, 5),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-6, 6))
      }else{
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature:",feature,"<br>","log2(FC): ",
                                    round(log2FC,2), "<br>", "-log10(padj) :",
                                    m_log10_p_adj), color=m_log10_p_adj,
                         size=2.5)) +
          ggplot2::guides(size="none") +
          ggplot2::labs(colour="-log10(padj)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0) +
          ggplot2::scale_y_continuous(breaks=c(-5, -2, -1, 0, 1, 2, 5),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-6, 6))
      }
      if(plotly == TRUE){
        if(nrow(sig.diff.exp) < 26){
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=500,
                                          tooltip="text")
        }else{
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=nrow(sig.diff.exp)*20,
                                          tooltip="text")
        }

        i <- NULL
        for(i in seq_len(length(in.dotchart$x$data))){
          in.dotchart$x$data[[i]]$text <- in.dotchart$x$data[[i]]$text %>%
            stringr::str_replace_all(pattern='log2(FC): 5',
                                     replacement='log2(FC): Inf') %>%
            stringr::str_replace_all(pattern='log2(FC):  5',
                                     replacement='log2(FC):  Inf') %>%
            stringr::str_replace_all(pattern='log2(FC): -5',
                                     replacement='log2(FC): -Inf') %>%
            stringr::str_replace_all("Inf\\..*", "Inf")
        }
      }else{
        in.dotchart <- p.dotchart
      }


    }else if(sum(is.infinite(sig.diff.exp$log2FC)) == 0){
      ##all sig. lipid is finite

      p.dotchart <- ggpubr::ggdotchart(sig.diff.exp, x="feature", y="log2FC",
                               rotate=TRUE,
                               color="white",
                               sorting="descending",
                               add="segments",
                               dot.size=2.5,
                               legend.title="log2(FC)",
                               xlab=" ",
                               ylab="log2(Fold Change)",
                               legend="right",
                               ggtheme=ggpubr::theme_pubr())
      p.dotchart <- if(sum(!is.na(sig.diff.exp[,'p_value'])) == 0){
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature, "<br>",
                                    "log2(FC) : ", round(log2FC, 2), "<br>",
                                    "-log10(p-value) :", m_log10_p_value),
                         color=log2FC,size=2.5))+
          ggplot2::guides(size='none') + ggplot2::labs(colour="log2(FC)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0)  +
          ggplot2::scale_y_continuous(breaks=c(-5, -2, -1, 0, 1, 2, 5),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-6, 6))
      }else if(sig_stat == "p"){
        if(nrow(sig.diff.exp) == 1){
          p.dotchart +
            ggplot2::geom_point(
              ggplot2::aes(text=paste("feature :", feature, "<br>",
                                      "log2(FC) : ", round(log2FC, 2), "<br>",
                                      "-log10(padj) :", m_log10_p_adj),
                           color=ifelse(m_log10_p_adj > 0, "red",
                                          "steelblue"),  dsize=2.5)) +
            ggplot2::guides(size="none") +
            ggplot2::theme(legend.position="none")
        }else{
          p.dotchart +
            ggplot2::geom_point(
              ggplot2::aes(text=paste("feature :", feature, "<br>",
                                      "log2(FC) : ", round(log2FC, 2), "<br>",
                                      "-log10(p-value) :", m_log10_p_value),
                           color=m_log10_p_value,size=2.5)) +
            ggplot2::guides(size="none") +
            ggplot2::labs(colour="-log10(p-value)") +
            ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                            high="red", midpoint=0)
        }
      }else{
        if(nrow(sig.diff.exp) == 1){
          p.dotchart +
            ggplot2::geom_point(
              ggplot2::aes(text=paste("feature :", feature, "<br>",
                                      "log2(FC) : ", round(log2FC, 2), "<br>",
                                      "-log10(padj) :", m_log10_p_adj),
                           color=ifelse(m_log10_p_adj > 0, "red",
                                        "steelblue"), size=2.5)) +
            ggplot2::guides(size="none") +
            ggplot2::theme(legend.position="none")
        }else{
          p.dotchart +
            ggplot2::geom_point(
              ggplot2::aes(text=paste("feature :", feature, "<br>",
                                      "log2(FC) : ", round(log2FC, 2), "<br>",
                                      "-log10(padj) :", m_log10_p_adj),
                           color=m_log10_p_adj,size=2.5)) +
            ggplot2::guides(size="none") +
            ggplot2::labs(colour="-log10(padj)") +
            ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                            high="red", midpoint=0)
        }
      }
      if(plotly == TRUE){
        if(nrow(sig.diff.exp) < 26){
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=500,
                                          tooltip="text")
        }else{
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=nrow(sig.diff.exp)*20,
                                          tooltip="text")
        }
      }else{
        in.dotchart <- p.dotchart
      }

    }else{

      rm.inf <- sig.diff.exp %>% dplyr::filter(is.finite(log2FC))
      INF <- ceiling(max(rm.inf$log2FC, na.rm=TRUE) * 1.5)

      sig.diff.exp$log2FC <- ifelse(is.infinite(sig.diff.exp$log2FC) &
                                      sig.diff.exp$log2FC > 0, INF,
                                    ifelse(is.infinite(sig.diff.exp$log2FC) &
                                             sig.diff.exp$log2FC < 0,
                                           -INF, sig.diff.exp$log2FC))

      p.dotchart <- ggpubr::ggdotchart(sig.diff.exp, x="feature", y="log2FC",
                               rotate=TRUE,
                               color="white",
                               sorting="descending",
                               add="segments",
                               dot.size=2.5,
                               legend.title="log2(FC)",
                               xlab=" ",
                               ylab="log2(Fold Change)",
                               legend="right",
                               ggtheme=ggpubr::theme_pubr()) #+
      p.dotchart <- if(sum(!is.na(sig.diff.exp[,'p_value'])) == 0){
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature, "<br>",
                                    "log2(FC) : ", round(log2FC, 2), "<br>",
                                    "-log10(p-value) :", m_log10_p_value),
                         color=log2FC, size=2.5))+
          ggplot2::guides(size='none') + ggplot2::labs(colour="log2(FC)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0)  +
          ggplot2::scale_y_continuous(breaks=c(-5, -2, -1, 0, 1, 2, 5),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-6, 6))
      }else if(sig_stat == "p"){
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature, "<br>",
                                    "log2(FC) : ", round(log2FC, 2),"<br>",
                                    "-log10(p-value) :", m_log10_p_value),
                         color=m_log10_p_value, size=2.5))+
          ggplot2::guides(size="none") +
          ggplot2::labs(colour="-log10(p-value)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0) +
          ggplot2::scale_y_continuous(breaks=c(-INF, -2, -1, 0, 1, 2, INF),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-(INF+1), (INF+1)))
      }else{
        p.dotchart +
          ggplot2::geom_point(
            ggplot2::aes(text=paste("feature :", feature, "<br>",
                                    "log2(FC) : ", round(log2FC,2), "<br>",
                                    "-log10(padj) :", m_log10_p_adj),
                         color=m_log10_p_adj, size=2.5)) +
          ggplot2::guides(size="none") +
          ggplot2::labs(colour="-log10(padj)") +
          ggplot2::scale_colour_gradient2(low="steelblue", mid="white",
                                          high="red", midpoint=0) +
          ggplot2::scale_y_continuous(breaks=c(-INF, -2, -1, 0, 1, 2, INF),
                                      labels=c('-Inf', -2, -1, 0, 1, 2, 'Inf'),
                                      limits=c(-(INF+1), (INF+1)))
      }
      if(plotly == TRUE){
        if(nrow(sig.diff.exp) < 26){
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=500,
                                          tooltip="text")
        }else{
          in.dotchart <- plotly::ggplotly(p.dotchart,
                                          height=nrow(sig.diff.exp)*20,
                                          tooltip="text")
        }
        i <- NULL
        for(i in seq_len(length(in.dotchart$x$data))){
          in.dotchart$x$data[[i]]$text <- in.dotchart$x$data[[i]]$text %>%
            stringr::str_replace_all(pattern=paste0('log2(FC): ', INF),
                                     replacement='log2(FC): Inf') %>%
            stringr::str_replace_all(pattern=paste0('log2(FC):  ', INF),
                                     replacement='log2(FC):  Inf') %>%
            stringr::str_replace_all(pattern=paste0('log2(FC): -', INF),
                                     replacement='log2(FC): -Inf') %>%
            stringr::str_replace_all("Inf\\..*", "Inf")
        }
      }else{
        in.dotchart <- p.dotchart
      }

    }

  }else{

    in.dotchart <- NULL

  } #if(nrow(sig.diff.exp) > 0)

  plot_MA_Vol <-if(sum(!is.na(exp_data_stat[,'p_value'])) == 0){
    exp_data_stat %>%
      dplyr::select(feature, mean_exp, mean_ctrl,
                    p_value, p_adj, 'M'='log2FC') %>%
      dplyr::mutate(A=(log2(mean_exp) + log2(mean_ctrl))/2,
             m.log.p=NA,
             sig_log2fc.pvalue_colors=if (sig_stat == "p"){
               ifelse( M > log2(sig_FC) ,'up-regulated',
                       ifelse(-M > log2(sig_FC),'down-regulated','none'))
             }else{
               ifelse( M > log2(sig_FC), 'up-regulated',
                       ifelse(-M > log2(sig_FC),'down-regulated','none'))
             })
  }else{
    exp_data_stat %>%
      dplyr::select(feature, mean_exp, mean_ctrl,
                    p_value, p_adj, 'M'='log2FC') %>%
      dplyr::mutate(A=(log2(mean_exp) + log2(mean_ctrl))/2,
             m.log.p=-log10(p_value),
             sig_log2fc.pvalue_colors=if(sig_stat == "p"){
               ifelse( M>log2(sig_FC) &
                         -log10(p_value) > -log10(sig_pvalue), 'up-regulated',
                       ifelse(-M>log2(sig_FC) &
                                -log10(p_value) > -log10(sig_pvalue),
                              'down-regulated', 'none'))
             }else{
               ifelse( M > log2(sig_FC) &
                         -log10(p_adj) > -log10(sig_pvalue), 'up-regulated',
                       ifelse(-M > log2(sig_FC) &
                                -log10(p_adj) > -log10(sig_pvalue),
                              'down-regulated', 'none'))
             })
  }
  plot_MA_Vol$sig_log2fc.pvalue_colors <-
    factor(plot_MA_Vol$sig_log2fc.pvalue_colors,
           levels=c("down-regulated", "none", "up-regulated"))
  if(plotly == TRUE){
    DE.species.maplot <- plotly::plot_ly(data=plot_MA_Vol,
                                         x=~as.numeric(A),
                                         y=~as.numeric(M),
                                         type="scatter",
                                         mode="markers",
                                         color=~ sig_log2fc.pvalue_colors,
                                         colors=c("#4169E1", "#DDDDDD",
                                                  "#FF4500"),
                                         showlegend=TRUE,
                                         marker=list(size=4),
                                         hoverinfo="text",
                                         text=~ paste(
                                           "</br>Lipid:",
                                           plot_MA_Vol$feature,
                                           "</br>A value:",
                                           round(as.numeric(A), 4),
                                           "</br>M value:",
                                           round(as.numeric(M), 4))) %>%
      plotly::layout(
        xaxis=list(
          title="A = (log<sub>2</sub>(exp)+log<sub>2</sub>(ctrl))/2"),
        yaxis=list(title="M = log<sub>2</sub>(exp)-log<sub>2</sub>(ctrl)"),
        title="MA Plot",
        legend=list(title=list(text="log2FC Significant"),
                    orientation='h',
                    xanchor="center",
                    x=0.5,
                    y=-0.18))
  }else{
    DE.species.maplot <- ggplot2::ggplot(plot_MA_Vol,
                    ggplot2::aes(x=as.numeric(A),
                                 y=as.numeric(M),
                                 color=sig_log2fc.pvalue_colors)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values=c("#4169E1", "#DDDDDD", "#FF4500")) +
      ggplot2::labs(color='log2FC Significant',
                    y='M = log2(exp)-log2(ctrl)',
                    x="A = (log2(exp)+log2(ctrl))/2",
                    title="MA Plot") +
      ggthemes::theme_hc() +
      ggplot2::geom_hline(yintercept=0) +
      ggplot2::geom_vline(xintercept=0)
  }

  if(sum(!is.na(exp_data_stat[,'p_value'])) == 0){
    DE.species.volcano <- NULL
  }else{
    if(plotly == TRUE){
      DE.species.volcano <- if (sig_stat == "p"){
        plotly::plot_ly(data=plot_MA_Vol,
                        x=~as.numeric(M),
                        y=~as.numeric(m.log.p),
                        type="scatter",
                        mode="markers",
                        color=~sig_log2fc.pvalue_colors,
                        colors=c("#4169E1", "#DDDDDD", "#FF4500"),
                        marker=list(size=4),
                        hoverinfo="text",
                        text=~ paste(
                          "</br>Lipid:",
                          plot_MA_Vol$feature,
                          "</br>A value:",
                          round(as.numeric(A), 4),
                          "</br>M value:",
                          round(as.numeric(M), 4),
                          "</br>-log10(p-value):",
                          round(-log10(p_value), 4),
                          "</br>-log10(padj):",
                          round(-log10(p_adj), 4))) %>%
          plotly::layout(
            xaxis=
              list(title=
                     "M = log<sub>2</sub>(exp)-log<sub>2</sub>(ctrl)"),
            yaxis=list(title="-log<sub>10</sub>(p-value)"),
            title="Volcano Plot",
            legend=list(title=list(text="Significant lipid"),
                        orientation='h',
                        xanchor="center",
                        x=0.5,
                        y=-0.18)
          )
      }else{plotly::plot_ly(data=plot_MA_Vol,
                            x=~as.numeric(M),
                            y=~as.numeric(-log10(p_adj)),
                            type="scatter",
                            mode="markers",
                            color=~sig_log2fc.pvalue_colors,
                            colors=c("#4169E1", "#DDDDDD", "#FF4500"),
                            marker=list(size=4),
                            hoverinfo="text",
                            text=~ paste(
                              "</br>Lipid:",
                              plot_MA_Vol$feature,
                              "</br>A value:",
                              round(as.numeric(A), 4),
                              "</br>M value:",
                              round(as.numeric(M), 4),
                              "</br>-log10(p-value):",
                              round(-log10(p_value), 4),
                              "</br>-log10(padj):",
                              round(-log10(p_adj), 4))) %>%
          plotly::layout(
            xaxis=list(title="M = log<sub>2</sub>(exp)-log<sub>2</sub>(ctrl)"),
            yaxis=list(title="-log<sub>10</sub>(padj)"),
            title="Volcano Plot",
            legend=list(title=list(text="Significant lipid"),
                        orientation='h',
                        xanchor="center",
                        x=0.5,
                        y=-0.18)
          )
      }
    }else{
      DE.species.volcano <- if (sig_stat == "p"){
        ggplot2::ggplot(plot_MA_Vol,
                        ggplot2::aes(x=as.numeric(M),
                                     y=as.numeric(m.log.p),
                                     color=sig_log2fc.pvalue_colors)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(values=c("#4169E1", "#DDDDDD",
                                               "#FF4500")) +
          ggplot2::labs(color='Significant lipid',
                        y='-log10(p-value)',
                        x="M = log2(exp)-log2(ctrl)",
                        title="Volcano Plot") +
          ggthemes::theme_hc() +
          ggplot2::geom_hline(yintercept=0) +
          ggplot2::geom_vline(xintercept=0)
      }else{
        ggplot2::ggplot(plot_MA_Vol,
                        ggplot2::aes(x=as.numeric(M),
                                     y=as.numeric(-log10(p_adj)),
                                     color=sig_log2fc.pvalue_colors)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(values=c("#4169E1", "#DDDDDD",
                                               "#FF4500")) +
          ggplot2::labs(color='Significant lipid',
                        y='-log10(padj)',
                        x="M = log2(exp)-log2(ctrl)",
                        title="Volcano Plot") +
          ggthemes::theme_hc() +
          ggplot2::geom_hline(yintercept=0) +
          ggplot2::geom_vline(xintercept=0)
      }
    }

  }

  return(list(DE_species_table_all=exp_data_stat,
              DE_species_table_sig=sig.diff.exp,
              DE_species_dotchart_sig=in.dotchart,
              DE_species_maplot=DE.species.maplot,
              DE_species_volcano=DE.species.volcano))

} #function


