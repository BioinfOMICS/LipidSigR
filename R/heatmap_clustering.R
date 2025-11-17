#' @title heatmap_clustering
#' @description Hierarchical clustering of lipid species or characteristics
#' derived from two/multiple groups.
#' @param de_se The resulting SummarizedExperiment object from the differential
#' expression analysis function, such as \code{\link{deSp_twoGroup}}, \code{\link{deSp_multiGroup}},
#' \code{\link{deChar_twoGroup}}, and \code{\link{deChar_multiGroup}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @param distfun Character. The distance measure for computing correlation coefficient
#' (or covariance). Allowed methods include "pearson", "kendall", and "spearman".
#' Default is \code{'pearson'}
#' @param hclustfun Character. The agglomeration method. This should be
#' (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete",
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC), or "centroid" (= UPGMC).
#' Default is \code{'complete'}
#' @param type Character. It must be 'all' or 'sig.' 'all' for output the results of all
#' lipid species or characteristics, and 'sig' for output the results of
#' significant lipid species or characteristics.
#' @return Return a list with 1 figure and 1 matrix.
#' \enumerate{
#' \item interactive_heatmap & static_heatmap: a heatmap provides an overview of
#' all/significant lipid species or characteristics that illustrates the differences between groups.
#' \item corr_coef_matrix: the matrix of the heatmap.
#' }
#' @export
#' @examples
#' data("de_data_twoGroup")
#' processed_se <- data_process(
#'     de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' deSp_se <- deSp_twoGroup(processed_se, ref_group='ctrl', test='t-test',
#'     significant='padj', p_cutoff=0.05, FC_cutoff=2, transform='log10')
#' result <- heatmap_clustering(de_se=deSp_se, char='class', distfun='pearson',
#'     hclustfun='complete', type='sig')

heatmap_clustering <- function(
        de_se, char, distfun=c('pearson', 'kendall', 'spearman'),
        hclustfun=c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'),
        type=c('all', 'sig') ){
    ## check data
    .check_de_outputSE(de_se, de_type="all")
    if(!is.null(S4Vectors::metadata(de_se)[["char"]])) {
        char <- S4Vectors::metadata(de_se)[["char"]]
        message("char ", char, " has been selected in upstream function.")
    } else {
        if (is.null(char) | isFALSE(.check_char(de_se, char, type='common'))) {
            stop("Wrong char input, you can view the available char list by list_lipid_char function.")
        }
    }
    if (is.null(distfun) | isFALSE(distfun %in% c('pearson', 'kendall', 'spearman')) ) {
        stop("distfun must be one of 'pearson', 'kendall', or 'spearman'.")
    }
    if (is.null(hclustfun) | isFALSE(hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')) ) {
        stop("hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    }
    if (is.null(type) | isFALSE(type %in% c('all', 'sig')) ) {
        stop("type must be one of 'all' or 'sig'.")
    }

    abundance_data <- .extract_df(de_se, type = "abundance")
    rownames(abundance_data) <- NULL
    group_info <- .extract_df(de_se, type = "group")
    rownames(group_info) <- NULL
    if(!is.null(S4Vectors::metadata(de_se)[["char"]])) {
        DE_result_table <- S4Vectors::metadata(de_se)[["sig_deChar_result"]]
        de_type <- "char"
    } else {
        DE_result_table <- S4Vectors::metadata(de_se)[["sig_deSp_result"]]
        de_type <- "sp"

    }
    lipid_char_table <- .extract_df(de_se, "lipid")
    if (nrow(abundance_data) <= 1)  {
        stop("Insufficient number of lipids.")
    }
    if (isTRUE(is.character(DE_result_table) && type == 'sig') ){
        stop("Insufficient number of significant lipids.")
    }
    #colnames(DE_result_table)[1] <- 'feature' ## to be check
    abundance.mat <- abundance_data %>%
        dplyr::select(feature, group_info$sample_name)
    if(type=='sig'){
        abundance.mat <- abundance.mat %>%
            dplyr::filter(feature %in% DE_result_table$feature)
    }
    if (nrow(abundance.mat) <= 2)  {
        stop("Insufficient number of lipids.")
    }
    abundance.mat %<>% tibble::column_to_rownames(var='feature') %>% as.matrix()
    colnames(abundance.mat) <- group_info$label_name
    abundance.mat <- sweep(abundance.mat, 1, rowMeans(abundance.mat, na.rm=TRUE))
    abundance.mat <- sweep(abundance.mat, 1, apply(abundance.mat, 1, sd, na.rm=TRUE), "/")
    if(sum(is.na(abundance.mat)) > 0){
        abundance.mat <- abundance.mat[-which(is.na(abundance.mat[, 2])), ]
    }
    ## create colGroup data frame
    colGroup <- data.frame(Sample=group_info$group, stringsAsFactors=FALSE)

    if (de_type=="char") {
        row_color_label <- NULL
        rowGroup <- NULL
    } else {
        rowGroup <- abundance.mat %>% as.data.frame() %>%
            tibble::rownames_to_column(var='feature') %>% dplyr::select(feature) %>%
            dplyr::left_join(lipid_char_table, by='feature') %>%
            dplyr::select(tidyselect::all_of(char)) %>%
            dplyr::pull()
        row_color <- grDevices::rainbow(length(unique(rowGroup)))
        row_color_label <- rowGroup
        for(j in seq(unique(rowGroup))){
            row_color_label[which(rowGroup == unique(rowGroup)[j])] <- row_color[j]
        }
        rowGroup[is.na(rowGroup)] <- ''
    }

    ## heatmap
    heat.res <- suppressWarnings(
        .heatmap_in(
            abundance.mat, lipid_char_table, char, distfun, hclustfun, colGroup, rowGroup, row_color_label)
    )
    static_heatmap <- suppressWarnings(
        .heatmap_static(
            abundance.mat, lipid_char_table, char, distfun, hclustfun, colGroup, rowGroup)
    )
    return(list(
        interactive_heatmap=heat.res$heatmap, static_heatmap=static_heatmap,
        corr_coef_matrix=heat.res$reorder.data))
}

.heatmap_in <- function(
        data, lipid_char_table, char, distfun, hclustfun,
        colGroup, rowGroup, row_color_label){
    all_col_text_size <- .col_text_size(data)
    all_row_text_size <- .row_text_size(data)
    data <- sweep(data, 1, rowMeans(data, na.rm=TRUE))
    data <- sweep(data, 1, apply(data, 1,sd, na.rm=TRUE), "/")
    if(sum(is.na(data)) > 0){
        data <- data[-which(is.na(data[, 2])), ]
    }
    cb_grid <- iheatmapr::setup_colorbar_grid(y_length=0.6, x_start=1, y_start=0.4)
    if(distfun %in% c("pearson", "kendall", "spearman")){
        col_dend <- stats::hclust(
            stats::as.dist(1-stats::cor(data, method=distfun)), method=hclustfun)
        row_dend <- stats::hclust(
            stats::as.dist(1-stats::cor(t(data), method=distfun)), method=hclustfun)
    }else{
        col_dend <- stats::hclust(stats::dist(t(data), method=distfun),method=hclustfun)
        row_dend <- stats::hclust(stats::dist(data, method=distfun), method=hclustfun)
    }
    if(min(data) >= 0 || max(data) <= 0){
        heatmap <- iheatmapr::iheatmap(data, colors=.heatmap_color_scale(data),
                                       colorbar_grid=cb_grid, scale="rows")
    }else{
        heatmap <- iheatmapr::iheatmap(data, colorbar_grid=cb_grid, scale="rows")
    }
    heatmap %<>% iheatmapr::add_col_annotation(
            annotation=colGroup, side="top", show_colorbar=FALSE)
    if(!is.null(char)){
        heatmap %<>% iheatmapr::add_row_annotation(
            annotation=rowGroup, colors=row_color_label, side="right",
            show_colorbar=FALSE)
    }
    heatmap %<>% iheatmapr::add_col_dendro(col_dend, side="top", reorder=TRUE) %>%
        iheatmapr::add_row_dendro(row_dend, side="right", reorder =TRUE)
    if(ncol(data) <= 50){
        heatmap %<>% iheatmapr::add_col_labels(side="bottom", size=all_col_text_size)
    }
    if(nrow(data)<=50){
        heatmap %<>% iheatmapr::add_row_labels(side="left", size=all_row_text_size)
    }
    reorder.data <- data[rev(row_dend$order), col_dend$order]
    return(list(heatmap=heatmap, reorder.data=reorder.data))
}

.heatmap_static <- function(
        abundance.mat, lipid_char_table, char, distfun, hclustfun,
        colGroup, rowGroup){
    if(!is.null(char) && !is.null(rowGroup)){
        rightAnn <- ComplexHeatmap::rowAnnotation(
            annotation=rowGroup,
            col=list(annotation=.colorAnn(rowGroup, type='row')),
            show_legend=FALSE)
    }else{
        rightAnn <- NULL
    }
    colName <- ifelse(ncol(abundance.mat) <= 50, TRUE, FALSE)
    rowName <- ifelse(nrow(abundance.mat) <= 50, TRUE, FALSE)
    heatmap <- ComplexHeatmap::Heatmap(
        matrix=abundance.mat, row_names_side="left", row_dend_side='right',
        column_names_side="bottom", right_annotation=rightAnn,
        col=.colorScale(abundance.mat), name='Signal',
        top_annotation=ComplexHeatmap::HeatmapAnnotation(
            sample=colGroup$Sample,
            col=list(sample=.colorAnn(colGroup$Sample, type='col')),
            show_legend=FALSE),
        show_row_names=rowName, show_column_names=colName,
        clustering_distance_rows=distfun, clustering_method_rows=hclustfun,
        clustering_distance_columns=distfun,
        clustering_method_columns=hclustfun)
    return(heatmap)
}

.colorAnn <- function(annotation,type=c('row', 'col')){
    annotation <- unique(annotation)
    nAnn <- length(annotation)
    if(type == 'row'){
        if(nAnn <= 2){
            colorPalette <- RColorBrewer::brewer.pal(n=nAnn, name="Paired")
            colorPalette <- colorPalette[seq(nAnn)]
        }else if(nAnn <= 12){
            colorPalette <- RColorBrewer::brewer.pal(n=nAnn, name="Paired")
        }else{
            colorPalette <- randomcoloR::distinctColorPalette(nAnn)
        }
        colorPalette <- stats::setNames(colorPalette, annotation)
    }else{
        if(nAnn <= 2){
            colorPalette <- RColorBrewer::brewer.pal(n=nAnn, name="Dark2")
            colorPalette <- colorPalette[seq(nAnn)]
        }else if(nAnn <= 9){
            colorPalette <- RColorBrewer::brewer.pal(n=nAnn, name="Dark2")
        }else{
            colorPalette <- randomcoloR::distinctColorPalette(nAnn)
        }
        colorPalette <- stats::setNames(colorPalette, annotation)
    }
    return(colorPalette)
}

.row_text_size <- function(data){
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
    return(all_row_text_size=all_row_text_size)
}
.col_text_size <- function(data){
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
    return(all_col_text_size=all_col_text_size)
}
