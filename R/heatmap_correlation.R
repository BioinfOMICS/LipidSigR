#' @title heatmap_correlation
#' @description The correlation heatmap illustrates the correlation between
#' lipid classes or samples and depicts the patterns in each group.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list
#' returned by \code{\link{list_lipid_char}}.
#' @param transform Character. Method for data transformation. Allowed methods
#' include "none", "log10", "square", and "cube". Select 'none' to skip data transformation.
#' Default is \code{'log10'}.
#' @param correlation Character. The method for computing correlation coefficient.
#' Allowed methods includes "pearson" and "spearman". Default is \code{'pearson'}.
#' @param distfun Character. The distance measure for computing correlation
#' coefficient (or covariance). Allowed methods include "pearson", "spearman",
#' "kendall", euclidean", "maximum", "manhattan", "canberra", "binary",
#' "minkowski". Default is \code{'maximum'}.
#' @param hclustfun Character. The agglomeration method. This should be
#' (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#' "median" (= WPGMC) or "centroid" (= UPGMC). Default is \code{'average'}.
#' @param type Character. It must be 'sample' or 'class'. 'sample' outputs the
#' correlation results of samples, and 'class' outputs output the correlation
#' results of lipid class.
#' @return Return a list of 1 interactive plot, 1 static plot, and 1 table.
#' \enumerate{
#' \item interactive_heatmap & static_heatmap: heatmaps illustrate the correlation between samples or lipid class
#' \item corr_coef_matrix: the matrix of the heatmap.
#' }
#' @export
#' @examples
#' data("profiling_data")
#' processed_se <- data_process(
#'     profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' result <- heatmap_correlation(processed_se, char=NULL, transform='log10',
#'     correlation='pearson', distfun='maximum', hclustfun='average', type='sample')

heatmap_correlation <- function(
        processed_se, char, transform=c('none', 'log10', 'square', 'cube'),
        correlation=c('pearson', 'spearman'),
        distfun=c('pearson', 'spearman', 'kendall', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'),
        hclustfun=c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'),
        type=c('sample', 'class')){

    ## check parameter
    .check_inputSE(processed_se, metadata_list=NULL)
    if (type=='class' && is.null(char)) {
        stop("For 'class' type, char cannot be NULL.")
    } else if (isFALSE(.check_char(processed_se, char, type='common'))) {
        stop("Wrong char input, you can view the available char list by list_lipid_char function.")
    }
    if (is.null(transform) | isFALSE(transform %in% c('none', 'log10', 'square', 'cube')) ) {
        stop("transform must be one of 'none', 'log10', 'square', or 'cube'.")
    }
    if (is.null(correlation) | isFALSE(correlation %in% c('pearson', 'spearman')) ) {
        stop("correlation must be one of 'pearson' or 'spearman'.")
    }
    if (is.null(distfun) | isFALSE(distfun %in% c('pearson', 'spearman', 'kendall', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski')) ) {
        stop("distfun must be one of 'pearson', 'spearman', 'kendall', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', or 'minkowski'.")
    }
    if (is.null(hclustfun) | isFALSE(hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')) ) {
        stop("hclustfun must be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid'.")
    }
    if (is.null(type) | isFALSE(type %in% c('sample', 'class')) ) {
        stop("type must be one of 'sample' or 'class'.")
    }

    abundance_raw <- .extract_df(processed_se, type="abundance")
    lipid_char_table <- .extract_df(processed_se,type='lipid')

    .check_imputation(abundance_raw)

    if(type == 'sample'){
        abundance <- .transform(abundance_raw, transform)
        sample_abundance <- subset(abundance, select=-feature)
        corr_coef <- stats::cor(
            sample_abundance, method=correlation, use="pairwise.complete.obs")
        ## add not showing warnings
        corr_p <- suppressWarnings(.cor.test.p(x=sample_abundance, method=correlation))
        if(all(!is.na(corr_coef)) & nrow(corr_coef)>=2 & ncol(corr_coef)>=2){
            res <- .sample_heatmap(corr_coef, sample_abundance, distfun, hclustfun)
            reorder_corr_coef <- res$reorder_corr_coef
            heatmap <- res$sample_hm
            static_plot <- .static_heatmap(corr_coef, correlation, distfun, hclustfun)
        } else {
            stop("Not enough data for plotting the heatmap.")
        }

    } else if(type == 'class'){
        ### sum characteristics
        abundance_ga <- abundance_raw %>% dplyr::left_join(
            lipid_char_table[c('feature', char)], by='feature') %>%
            dplyr::select(-1)
        abundance_ga <- abundance_ga[!is.na(abundance_ga[[char]]),]
        if(nrow(abundance_ga) == 0){
            stop("Not enough data for plotting the heatmap.")
        } else {
            ## process "|" characteristic
            if(any(stringr::str_detect(abundance_ga[, char], '\\|'))  ){
                abundance_sum <- tidyr::separate_rows(abundance_ga, dplyr::all_of(char), sep ="\\|")
            } else {
                abundance_sum <- abundance_ga
            }
            abundance_sum <- abundance_sum %>%
                stats::aggregate(stats::as.formula(stringr::str_c('. ~ ',char)), ., sum)
        }

        abundance <- .transform(abundance_sum, transform)
        lipids_abundance <- abundance %>%
            tidyr::gather(-dplyr::all_of(char), key='sample_name', value='value') %>%
            tidyr::spread(key=char, value='value')
        lipids_abundance <- subset(lipids_abundance, select=-sample_name)
        corr_coef <- stats::cor(lipids_abundance, method=correlation,
                                use="pairwise.complete.obs")
        corr_p <- .cor.test.p(x=lipids_abundance, method=correlation)
        if(all(!is.na(lipids_abundance)) & nrow(lipids_abundance)>=2 & ncol(lipids_abundance)>=2){
            res <- .lipid_heatmap(corr_coef, lipids_abundance, distfun, hclustfun)
            reorder_corr_coef <- res$reorder_corr_coef
            heatmap <- res$lipids_hm
            static_plot <- .static_heatmap(corr_coef, correlation, distfun, hclustfun)
        } else {
            stop("Not enough data for plotting the heatmap.")
        }
    }
    return(list(
        interactive_heatmap=heatmap, static_heatmap=static_plot,
        corr_coef_matrix=reorder_corr_coef))
}

.cor.test.p <- function(x, method){
    FUN <- function(x, y, method) stats::cor.test(
        x=x, y=y, method=method, na.action="na.exclude")[["p.value"]]
    z <- outer(colnames(x), colnames(x),
               Vectorize(function(i, j) FUN(x[, i], x[, j], method=method)))
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}

.sample_heatmap <- function(corr_coef, sample_abundance, distfun, hclustfun){
    cb_grid <- iheatmapr::setup_colorbar_grid(
        y_length =0.6,x_start=1,y_start=0.5)
    if(distfun%in%c("pearson","kendall","spearman")){
        col_dend <- hclust(as.dist(
            1-cor(sample_abundance, method=distfun)), method=hclustfun)
    }else{
        col_dend <- hclust(dist(
            t(sample_abundance), method=distfun),method=hclustfun)
    }
    if(min(corr_coef)>=0 || max(corr_coef)<=0){
        sample_hm <- iheatmapr::iheatmap(
            corr_coef, colors=.heatmap_color_scale(corr_coef), colorbar_grid=cb_grid) %>%
            iheatmapr::add_col_dendro(
                col_dend,side="top",reorder =TRUE,size=0.1) %>%
            iheatmapr::add_row_dendro(
                col_dend,side="right",reorder =TRUE,size=0.1)
    }else{
        sample_hm <- iheatmapr::iheatmap(
            corr_coef,colorbar_grid=cb_grid) %>%
            iheatmapr::add_col_dendro(
                col_dend,side="top",reorder =TRUE,size=0.1) %>%
            iheatmapr::add_row_dendro(
                col_dend,side="right",reorder =TRUE,size=0.1)
    }
    if(nrow(corr_coef)<50){
        sample_hm <- sample_hm %>%
            iheatmapr::add_row_labels(font=list(size=10))
    }
    if(ncol(corr_coef)<50){
        sample_hm <- sample_hm %>%
            iheatmapr::add_col_labels(side="bottom",font=list(size=10))
    }
    ## reorder Sample correlation matrix
    reorder_corr_coef <- apply(corr_coef[,col_dend$order],2,rev)
    return(list(sample_hm=sample_hm, reorder_corr_coef=reorder_corr_coef))
}

.lipid_heatmap <- function(corr_coef, lipids_abundance, distfun, hclustfun) {
    cb_grid <- iheatmapr::setup_colorbar_grid(
        y_length =0.6,x_start=1,y_start=0.5)
    if(distfun%in%c("pearson","kendall","spearman")){
        col_dend <- hclust(
            as.dist(1-cor(lipids_abundance, method=distfun)),
            method=hclustfun)
    }else{
        col_dend <- hclust(
            dist(t(lipids_abundance), method=distfun),method=hclustfun)
    }
    if(min(corr_coef)>0 | max(corr_coef)<0){
        lipids_hm <- iheatmapr::iheatmap(
            corr_coef,colors=.heatmap_color_scale(corr_coef),
            colorbar_grid=cb_grid) %>%
            iheatmapr::add_col_dendro(
                col_dend,side="top", reorder=TRUE, size=0.1) %>%
            iheatmapr::add_row_dendro(
                col_dend,side="right", reorder=TRUE, size=0.1)
    }else{
        lipids_hm <- iheatmapr::iheatmap(
            corr_coef,colorbar_grid=cb_grid) %>%
            iheatmapr::add_col_dendro(
                col_dend,side="top", reorder=TRUE, size=0.1) %>%
            iheatmapr::add_row_dendro(
                col_dend,side="right", reorder=TRUE, size=0.1)
    }
    if(nrow(corr_coef)<50){
        lipids_hm <- lipids_hm %>%
            iheatmapr::add_row_labels(font=list(size=10))
    }
    if(ncol(corr_coef)<50){
        lipids_hm <- lipids_hm %>%
            iheatmapr::add_col_labels(side="bottom",font=list(size=10))
    }
    ## reorder Lipids correlation matrix
    reorder_corr_coef <- apply(corr_coef[,col_dend$order],2,rev)
    return(list(lipids_hm=lipids_hm, reorder_corr_coef=reorder_corr_coef))
}

.static_heatmap <- function(corr_coef, correlation, distfun, hclustfun){

    col_dend <- if(distfun %in% c("pearson", "kendall", "spearman")) {
        stats::hclust(
            stats::as.dist(
                1 - stats::cor(corr_coef, method=distfun)), method=hclustfun)
    } else {
        stats::hclust(
            stats::dist(t(corr_coef), method=distfun), method=hclustfun)
    }
    row_dend <- col_dend %>% stats::as.dendrogram() %>% rev()
    heatmap <- ComplexHeatmap::Heatmap(
        matrix=corr_coef, row_names_side="left", row_dend_side='right',
        column_names_side="bottom",  col=.colorScale(corr_coef),
        name='Signal', cluster_rows=row_dend,
        cluster_columns=col_dend)
    return(heatmap)
}
