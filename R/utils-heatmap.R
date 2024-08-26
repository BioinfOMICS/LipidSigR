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
