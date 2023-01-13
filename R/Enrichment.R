#' @title Enrichment
#' @description
#' \enumerate{
#' \item Determine whether significant lipid species are enriched in
#' user-selected characteristics.
#' \item Print a summary table and enrichment bar chart.
#' }
#' @param DE_species_table_sig A data frame comprises the significant results
#' of differential expression analysis, including fold change, p-value,
#' adjusted p-value. The output of \code{\link{DE_species_2}}.
#' @param lipid_char_table A data frame with lipid features, such as class,
#' total length. NAs are allowed. The name of first column must be "feature".
#' @param char_var  A character string of the first lipid characteristic
#' selected by users from the column name of \bold{lipid_char_table},
#' such as total length.
#' @param sig_pvalue Numeric. Significant level. (default: 0.05)
#' @param plotly Logical value. If TRUE, return the resulting plots dynamically.
#' @importFrom magrittr %>%
#' @return Return a list with 1 tibble and 1 plot.
#' \enumerate{
#' \item enrich_char_table: a data frame encompasses condition, characteristic,
#'  p-value, and significance.
#' \item enrich_char_barplot: the enrichment plot for
#' up/down/non-significant groups
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
#'                                       replace_zero=TRUE, zero2what='min',
#'                                       xmin=0.5, replace_NA=TRUE,
#'                                       NA2what='min', ymin=0.5,
#'                                       pct_transform=TRUE,
#'                                       data_transform=FALSE,
#'                                       trans_type='log', centering=FALSE,
#'                                       scaling=FALSE)
#' lipid_char_filter <- lipid_char_table %>%
#'     filter(feature %in% exp_transform$feature)
#' DE_species_table_sig <- DE_species_2(exp_transform_non_log,
#'                                      data_transform=TRUE,
#'                                      group_info=group_info, paired=FALSE,
#'                                      test='t.test', adjust_p_method='BH',
#'                                      sig_stat='p.adj', sig_pvalue=0.05,
#'                                      sig_FC=2,
#'                                      plotly=TRUE)$DE_species_table_sig
#' char_var <- colnames(lipid_char_filter)[-1]
#' Enrichment (DE_species_table_sig, lipid_char_table=lipid_char_filter,
#'             char_var=char_var[1], sig_pvalue=0.05, plotly=TRUE)
Enrichment <- function(DE_species_table_sig, lipid_char_table, char_var,
                       sig_pvalue=0.05, plotly=TRUE){

  DE_species_table_sig <- as.data.frame(DE_species_table_sig)
  lipid_char_table <- as.data.frame(lipid_char_table)
  if(nrow(DE_species_table_sig) < 5){
    stop("Less than 5 significant lipid features")
  }
  if(!is(lipid_char_table[, 1], 'character')){
    stop("lipid_char_table first column must contain a list of
         lipids names (features).")
  }
  if(tibble::is_tibble(lipid_char_table)){
    if(nrow(lipid_char_table) != nrow(unique(lipid_char_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }else{
    if(nrow(lipid_char_table) != length(unique(lipid_char_table[, 1]))){
      stop("The lipids name (features) must be unique")
    }
  }
  if("class" %in% colnames(lipid_char_table)){
    if(!is(lipid_char_table[, 'class'], 'character')){
      stop("lipid_char_table content of column 'class' must be characters")
    }
  }
  if("totallength" %in% colnames(lipid_char_table)){
    if(!class(lipid_char_table[, 'totallength']) %in% c("integer", "numeric")){
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
  if(!char_var %in% colnames(lipid_char_table)){
    stop("char_var must be included in the lipid_char_table.")
  }
  if(ncol(dplyr::select(lipid_char_table,
                        tidyselect::starts_with("FA_"))) == 0){
    warning("(OPTIONAL) lipid_char_table does not contain column
            names starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>%
      dplyr::select(feature, tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_", colnames(FA_lipid_char_table), value=TRUE)
    max_comma <- 0
    for(i in seq_len(length(FA_col))){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[, col],
                                            ','),na.rm=TRUE)
      if(comma_count > 0){
        FA_lipid_char_table <- FA_lipid_char_table %>%
          tidyr::separate(col,
                          c(col,
                            paste0(col, "_", seq_len(comma_count))),
                          ",",
                          convert=TRUE)
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>%
      tidyr::gather(lipid.category, lipid.category.value, -feature)
    if(max_comma > 0){
      for (i in seq_len(max_comma)) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <-FA_lipid_char_table[-intersect(
          grep(select_name,
               FA_lipid_char_table[, "lipid.category"]),
          which(is.na(FA_lipid_char_table$lipid.category.value))), ]
      }
    }
    if(is(FA_lipid_char_table$lipid.category.value, 'character')){
      stop("In the 'FA_' related analyses, the values are positive integer
           or zero and separated by comma. i.e., 10,12,11")
    }else if(sum(
      stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)) !=
      round(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value)))) != 0 |
      min(stats::na.omit(
        as.numeric(FA_lipid_char_table$lipid.category.value))) < 0){
      stop("In the 'FA_' related analyses, the values are positive integer
           or zero and separated by comma. i.e., 10,12,11")
    }
  }
  lipo.fisher <- function(sig.mat, lipo.category){

    sig.class <- sig.mat %>%
      dplyr::group_by(characteristic)%>%
      dplyr::summarise(sig.count=dplyr::n())%>%
      dplyr::left_join(lipo.category, by='characteristic')%>%
      dplyr::mutate(p.value=NA, significant='NO')

    n.sig <- sig.class%>%
      dplyr::summarise(sum=sum(sig.count)) %>%
      .$sum

    n.lipid <- lipo.category %>%
      dplyr::summarise(sum=sum(total.count)) %>%
      .$sum

    for (i in seq_len(nrow(sig.class))){

      sig.class$p.value[i] <- stats::fisher.test(
        matrix(c(sig.class$sig.count[i],
                 sig.class$total.count[i],
                 n.sig-sig.class$sig.count[i],
                 n.lipid-sig.class$total.count[i]), nrow=2),
        alternative='greater')$p.value
    }
    sig.class$m.log.p <- -log10(sig.class$p.value)
    sig.class$significant[sig.class$p.value < sig_pvalue] <- 'YES'

    return(sig.class)
  } #function


  lipo.category <- lipid_char_table %>%
    dplyr::select(feature, 'characteristic'=tidyselect::all_of(char_var)) %>%
    dplyr::group_by(characteristic) %>%
    dplyr::summarise(total.count=dplyr::n())

  char.tab <- lipid_char_table %>%
    dplyr::select(feature, 'characteristic'=tidyselect::all_of(char_var))

  sig_lipid.all <- DE_species_table_sig %>%
   dplyr::left_join(char.tab, by='feature') %>%
   dplyr::filter(!is.na(characteristic)) %>%
   dplyr::mutate(condition='UP & DOWN')

  sig_lipid.up <- DE_species_table_sig %>%
    dplyr::left_join(char.tab, by='feature') %>%
    dplyr::filter(!is.na(characteristic)) %>%
    dplyr::filter(log2FC > 0) %>%
    dplyr::mutate(condition='UP')

  sig_lipid.down <- DE_species_table_sig %>%
    dplyr::left_join(char.tab, by='feature') %>%
    dplyr::filter(!is.na(characteristic)) %>%
    dplyr::filter(log2FC < 0) %>%
    dplyr::mutate(condition='DOWN')

  #### enrichment analysis by category ####
  sig_class.all <- lipo.fisher(sig_lipid.all, lipo.category) %>%
    dplyr::mutate(condition='UP & DOWN')
  if(nrow(sig_lipid.up) > 0){
    sig_class.up <- lipo.fisher(sig_lipid.up, lipo.category) %>%
      dplyr::mutate(condition='UP')
  }else{
    sig_class.up <- NULL
  }
  if(nrow(sig_lipid.down) > 0){
    sig_class.down <- lipo.fisher(sig_lipid.down, lipo.category) %>%
      dplyr::mutate(condition='DOWN')
  }else{
    sig_class.down <- NULL
  }

  sig_class <- dplyr::bind_rows(sig_class.all,
                                sig_class.up, sig_class.down)
  sig_class$characteristic <- factor(sig_class$characteristic,
                                     sort(unique(sig_class$characteristic)))

  plot.tab <- sig_class %>%
    dplyr::filter(condition != 'UP & DOWN') %>%
    dplyr::mutate(mlogP=ifelse(condition == 'DOWN', -(m.log.p), m.log.p),
           significance=ifelse(significant == 'YES' & condition == 'UP', 'UP',
                                 ifelse(significant == 'YES' &
                                          condition == 'DOWN',
                                        'DOWN', 'non-significant')))

  if(max(plot.tab$m.log.p) == 0){

    in.sig.class <- NULL

  }else{

    if(length(unique(plot.tab$condition)) == 2){

      x.max <- max(ceiling(plot.tab$m.log.p))
      x.label <- as.character(c(seq(x.max,1), 0, seq(1, x.max)))
      plot.tab <- plot.tab %>%
        dplyr::group_by(characteristic) %>%
        dplyr::mutate(rank=mlogP[which(abs(mlogP) == max(abs(mlogP)))])
      p.sig.class <- ggplot2::ggplot(plot.tab,
                                     ggplot2::aes(x=mlogP,
                                                  y=stats::reorder(
                                                    characteristic,
                                                    rank, max),
                                                  fill=significance)) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept=0, color='#444444') +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(breaks=c('DOWN', 'UP', 'non-significant'),
                                   values=c('#4169E1', '#FF6347', '#666666'))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=14),
              axis.text.y=ggplot2::element_text(size=12),
              axis.text.x=ggplot2::element_text(size=12)) +
        ggplot2::labs(x='-log10(p-value)', y=char_var) +
        ggplot2::scale_x_continuous(breaks=-x.max:x.max,
                                    labels=x.label)

    }else if(length(unique(plot.tab$condition)) == 1 &
             unique(plot.tab$condition) == 'UP'){

      x.max <- max(ceiling(plot.tab$m.log.p))
      x.label <- as.character(0:x.max)

      p.sig.class <- ggplot2::ggplot(plot.tab,
                                     ggplot2::aes(x=mlogP,
                                                  y=stats::reorder(
                                                    characteristic,
                                                    mlogP),
                                                  fill=significance)) +
        ggplot2::geom_col() +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(breaks=c('DOWN', 'UP', 'non-significant'),
                                   values=c('#4169E1', '#FF6347', '#666666'))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=14),
              axis.text.y=ggplot2::element_text(size=12),
              axis.text.x=ggplot2::element_text(size=12)) +
        ggplot2::labs(x='-log10(p-value)', y=char_var) +
        ggplot2::scale_x_continuous(breaks=0:x.max,
                                    labels=x.label)

    }else if(length(unique(plot.tab$condition)) == 1 &
             unique(plot.tab$condition) == 'DOWN'){

      x.max <- max(ceiling(plot.tab$m.log.p))
      x.label <- as.character(x.max:0)

      p.sig.class <- ggplot2::ggplot(plot.tab,
                                     ggplot2::aes(x=mlogP,
                                                  y=stats::reorder(
                                                    characteristic,
                                                    mlogP),
                                                  fill=significance)) +
        ggplot2::geom_col() +
        ggthemes::theme_hc() +
        ggplot2::scale_fill_manual(breaks=c('DOWN', 'UP', 'non-significant'),
                                   values=c('#4169E1', '#FF6347', '#666666'))+
        ggplot2::theme(axis.title=ggplot2::element_text(size=14),
              axis.text.y=ggplot2::element_text(size=12),
              axis.text.x=ggplot2::element_text(size=12)) +
        ggplot2::labs(x='-log10(p-value)', y=char_var) +
        ggplot2::scale_x_continuous(breaks=-x.max:0,
                                    labels=x.label)

    }

    if(plotly == TRUE){
      p.sig.class$data$p.value <- round(p.sig.class$data$p.value, 3)
      p.sig.class$data$mlogP <- round(p.sig.class$data$mlogP, 3)
      in.sig.class <- plotly::ggplotly(p.sig.class)
      for (i in seq_len(length(in.sig.class$x$data))){
        in.sig.class$x$data[[i]]$text <-
          gsub("reorder\\(characteristic, rank, max\\)",
               char_var, in.sig.class$x$data[[i]]$text)
        in.sig.class$x$data[[i]]$text <-
          gsub("reorder\\(characteristic, mlogP\\)",
               char_var, in.sig.class$x$data[[i]]$text)
        in.sig.class$x$data[[i]]$text <-
          stringr::str_replace(string=in.sig.class$x$data[[i]]$text,
                               pattern='-', replacement='')
        in.sig.class$x$data[[i]]$text <-
          stringr::str_replace(string=in.sig.class$x$data[[i]]$text,
                               pattern='mlogP', replacement='-log10(p-value)')
        if (!is.null(in.sig.class$x$data[[i]]$hovertext)){
          in.sig.class$x$data[[i]]$hovertext <-
            gsub("\\(", "", in.sig.class$x$data[[i]]$hovertext)
          in.sig.class$x$data[[i]]$hovertext <-
            gsub("\\)", "", in.sig.class$x$data[[i]]$hovertext)
          in.sig.class$x$data[[i]]$hovertext <-
            gsub("as.charactersig.count", "sig.count ",
                 in.sig.class$x$data[[i]]$hovertext)
          in.sig.class$x$data[[i]]$hovertext <-
            gsub("\\reorder(characteristic, mlogP)",
                 char_var, in.sig.class$x$data[[i]]$hovertext)
          in.sig.class$x$data[[i]]$hovertext <-
            stringr::str_replace(in.sig.class$x$data[[i]]$hovertext,
                                 pattern='-', replacement='')
          in.sig.class$x$data[[i]]$hovertext <-
            stringr::str_replace(in.sig.class$x$data[[i]]$hovertext,
                                 pattern='mlogP',
                                 replacement='-log10(p-value)')
        }
      }
    }else{
      in.sig.class <- p.sig.class
    }


  }


  final.tab <- sig_class %>%
    dplyr::filter(condition != 'UP & DOWN') %>%
    dplyr::select(condition, characteristic,
                  sig.count, total.count,
                  p.value, m.log.p, significant)

  final.tab$m.log.p <- round(final.tab$m.log.p, 3)

  return(list(enrich_char_table=final.tab,
              enrich_char_barplot=in.sig.class))

} #function


