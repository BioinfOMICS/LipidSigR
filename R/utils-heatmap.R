.heatmap_color_scale <- function(data){
    data <- round(data, 3)
    if(max(data) <= 0 & min(data) < 0){
        over_median <- min(data)/2
        if(max(data) < over_median){
            color <-  grDevices::colorRampPalette(c("#157AB5", "#92c5de"))(n=1000)
        }else{
            color_rank <- round(max(data)/(min(data))*1000)
            color_scale <- grDevices::colorRampPalette(
                c("#0571b0", "#92c5de", "white"))(n = 1000)
            color <- color_scale[color_rank:1000]
        }
    }else if(min(data) >= 0 & max(data) > 0){
        over_median <- max(data)/2
        if(min(data) > over_median){
            color <-  grDevices::colorRampPalette(c("#f4a582", "#ca0020"))(n=1000)
        }else{
            color_rank <- round(min(data)/(max(data))*1000)
            color_scale <- grDevices::colorRampPalette(
                c("white", "#f4a582", "#ca0020"))(n = 1000)
            color <- color_scale[color_rank:1000]
        }
    }
    return(color)
}

.corr_heatmap <- function(
        table_all, table_sig, heatmap_col, lipid_char_table, side_color_char, significant,
        distfun, hclustfun, cor_type) {
    P.mat <- table_all %>%
        dplyr::select(clin_factor, feature, !!rlang::sym(significant)) %>%
        tidyr::spread(feature, !!rlang::sym(significant)) %>%
        tibble::column_to_rownames(var='clin_factor') %>% as.matrix()
    if(length(unique(table_sig[["feature"]])) > 1){
        max_colcex <- max(stringr::str_length(table_sig[["feature"]]))
        max_rowcex <- max(stringr::str_length(table_sig[["clin_factor"]]))
        if (max_colcex<4) {max_colcex <- 4}
        if (max_rowcex<4) {max_rowcex <- 4}

        Cor.mat <- table_sig %>%
            dplyr::select(clin_factor, feature, !!rlang::sym(heatmap_col)) %>%
            tidyr::spread(feature, !!rlang::sym(heatmap_col)) %>%
            tibble::column_to_rownames(var = 'clin_factor') %>% as.matrix()
        Cor.mat[is.na(Cor.mat)] <- 0

        if(sum(is.na(Cor.mat))==0 & nrow(Cor.mat) >= 2 & ncol(Cor.mat) >= 2){
            if (distfun %in% c("pearson","kendall","spearman")){
                col_dend <- hclust(as.dist(1-stats::cor(Cor.mat, method=distfun)), method=hclustfun)
                row_dend <- hclust(as.dist(1-stats::cor(t(Cor.mat), method=distfun)), method=hclustfun)
            }else{
                col_dend <- hclust(dist(t(Cor.mat), method=distfun), method=hclustfun)
                row_dend <- hclust(dist(Cor.mat, method=distfun), method=hclustfun)
            }
            ax <- list(title="", zeroline=FALSE, showline=FALSE, showgrid=FALSE, ticks='')
            bx <- list(title="", zeroline=FALSE, showline=FALSE, showticklabels=FALSE, showgrid=FALSE, ticks='')

            if (distfun %in% c('pearson','spearman','kendall')){
                dist_fun <- function(x){
                    x <- t(x)
                    if (cor_type == "cor") {
                        cor.mat <- stats::cor(x,method=distfun, use='complete.obs')
                    } else {
                        cor.mat <- stats::cor(x,method=distfun)
                    }
                    cor.mat <- (1-cor.mat)
                    cor.dist <- stats::as.dist(cor.mat)
                    return(cor.dist)
                }
            }else{
                dist_fun <- function(x) stats::dist(x, method=distfun)
            }

            hclust_fun <- function(x) stats::hclust(x, method=hclustfun)

            if (!is.null(side_color_char)){
                lipid_char_table <- lipid_char_table[match(colnames(Cor.mat), lipid_char_table$feature),]
                colGroup <- lipid_char_table %>% dplyr::select(all_of(side_color_char))
                hm <- Cor.mat %>% heatmaply::heatmapr(Rowv=row_dend, Colv=col_dend, col_side_colors=colGroup)
            }else{
                hm <- Cor.mat %>% heatmaply::heatmapr(Rowv=row_dend, Colv=col_dend)
            }
            in.heatmap <- hm %>% heatmaply::heatmaply(
                scale_fill_gradient_fun=ggplot2::scale_fill_gradient2(
                    low="#0571b0", mid="white", high="#ca0020", midpoint=0),
                node_type="heatmap", column_text_angle=270, grid_color="black",
                margins=c(l=0.2, r=0.8, t=20, b=80))

            if(nrow(Cor.mat)>50 & ncol(Cor.mat)>50){
                in.heatmap <- in.heatmap %>% plotly::layout(xaxis=bx, yaxis=bx)
            }else if(nrow(Cor.mat)>50){
                in.heatmap <- in.heatmap %>% plotly::layout(xaxis=ax, yaxis=bx)
            }else if(ncol(Cor.mat)>50){
                in.heatmap <- in.heatmap %>% plotly::layout(xaxis=bx, yaxis=ax)
            }else{
                in.heatmap <- in.heatmap %>% plotly::layout(xaxis=ax, yaxis=ax)
            }
            reorder_Cor.mat <- Cor.mat[row_dend$order,rev(col_dend$order)]
            P.mat <- P.mat[rev(rownames(reorder_Cor.mat)), colnames(reorder_Cor.mat)]
            for (i in seq(length(in.heatmap$x$data))) {
                if(!is.null(in.heatmap$x$data[[i]])){
                    if(length(in.heatmap$x$data[[i]]$text)==length(P.mat)){
                        for (j in seq(nrow(in.heatmap$x$data[[i]]$text))){
                            for(z in seq(ncol(in.heatmap$x$data[[i]]$text))){
                                in.heatmap$x$data[[i]]$text[j,z] <- paste0(
                                    in.heatmap$x$data[[i]]$text[j,z],'<br>P-value: ',P.mat[j,z])
                            }
                        }
                    }
                }
            }
            suppressWarnings(stats::heatmap(
                Cor.mat, Rowv=TRUE, Colv=TRUE, dendrogram='both', trace="none",
                col=grDevices::colorRampPalette(c("blue", "white", "red"))(n=300),
                distfun=dist_fun, hclustfun=hclust_fun, main=NULL, margins=c(8,8),
                key.xlab=heatmap_col, lwid=c(1, 9), cexCol=3/log2(max_colcex),
                cexRow=6/log2(max_rowcex), scale='none') )
            static_hm <- grDevices::recordPlot()
            grDevices::dev.off()
            reorder_Cor.mat <- Cor.mat[rev(row_dend$order), col_dend$order]
        } else{
            stop("Insufficient data for plotting heatmap.")
        }
    } else{
        stop("Insufficient data for plotting heatmap.")
    }
    return(list(in.heatmap=in.heatmap, static_hm=static_hm, plot.mat=reorder_Cor.mat))
}
