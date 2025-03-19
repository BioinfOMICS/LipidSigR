#' @title plot_data_process
#' @description This function plots the differences in abundance data before and after data processing.
#' @param se A SummarizedExperiment object construct by \code{\link{as_summarized_experiment}}.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @return Return a list with 4 static plots and 4 interactive plots.
#' \enumerate{
#' \item static_boxPlot_before: a static bar plot of lipid abundance before data processing.
#' \item static_densityPlot_before: a static density plot of lipid abundance before data processing.
#' \item static_boxPlot_after: a static bar plot of lipid abundance before data processing.
#' \item static_densityPlot_after: a static density plot of lipid abundance after data processing.
#' \item interactive_boxPlot_before: an interactive bar plot of lipid abundance before data processing.
#' \item interactive_densityPlot_before: an interactive density plot of lipid abundance before data processing.
#' \item interactive_boxPlot_after: an interactive bar plot of lipid abundance after data processing.
#' \item interactive_densityPlot_after: an interactive density plot of lipid abundance after data processing.
#' }
#' @export
#' @examples
#' data("profiling_data")
#' processed_se <- data_process(
#'     profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage', transform='log10')
#' plots <- plot_data_process(profiling_data, processed_se)
plot_data_process <- function(se, processed_se){
    ## check data

    ##
    abundance_before <- .extract_df(se, type="abundance")
    #abundance_after <- .extract_df(processed_se, type="abundance")
    abundance_after <- S4Vectors::metadata(processed_se)$processed_abund
    ## quality plot
    status_type <- c("before", "after")
    for (i in seq(status_type)) {
        data <- get0(paste0("abundance_", status_type[i]))
        status <- status_type[i]
        abundance_trans_data <- reshape2::melt(data, id.vars = colnames(data[1]))
        abundance_long <- stats::na.omit(abundance_trans_data)
        colnames(abundance_long) <- c("feature", "sample", "value")
        sample_list <- colnames(data)[-1][1:length(colnames(data)[-1])]
        sample_list_box <- sample_list[length(sample_list):1]
        ## box plot
        abundance_box <- abundance_long
        abundance_box$sample <- factor(abundance_box$sample, levels=sample_list_box)
        box <- ggplot2::ggplot(abundance_box, ggplot2::aes(x = sample , y = value)) +
            ggplot2::geom_boxplot(notch = TRUE, color="#09567A", fill="#c8dfea") +
            ggplot2::coord_flip() +
            ggplot2::theme_classic() +
            ggplot2::ylab("Abundance") + ggplot2::xlab("Samples") +
            ggplot2::scale_y_continuous(trans='log10') +
            ggplot2::labs(title = paste0("Box plot ", status, " data process"), subtitle = "") #+
        #stat_summary(fun = mean, geom = "point", size = 2, color = "#B78428")
        boxPlot <- plotly::ggplotly(box)
        ## density plot
        densityPlot <- .summaryDensity(abundance_long, sample_list, status)
        assign(paste0("interactive_boxPlot_", status), boxPlot)
        assign(paste0("static_boxPlot_", status), box)
        assign(paste0("interactive_densityPlot_", status), densityPlot$densityPlot)
        assign(paste0("static_densityPlot_", status), densityPlot$density)
    }

    return(list(
        interactive_boxPlot_before=interactive_boxPlot_before,
        static_boxPlot_before=static_boxPlot_before,
        interactive_densityPlot_before=interactive_densityPlot_before,
        static_densityPlot_before=static_densityPlot_before,
        interactive_boxPlot_after=interactive_boxPlot_after,
        static_boxPlot_after=static_boxPlot_after,
        interactive_densityPlot_after=interactive_densityPlot_after,
        static_densityPlot_after=static_densityPlot_after))
}

.summaryDensity <- function(abundance_long, sample_list, status){
    ## static
    abundance_dis <- abundance_long
    abundance_dis$sample <- factor(abundance_dis$sample, levels = sample_list)
    density <- ggplot2::ggplot(abundance_dis, ggplot2::aes(x=value, color=sample)) +
        ggplot2::geom_density() +
        ggplot2::theme_classic() +
        ggplot2::ylab("Density") + ggplot2::xlab("Abundance") +
        #ggplot2::scale_x_continuous(trans='log10') +
        #ylab("Density") + xlab("log(expression)") + xlim(c(-1,4)) +
        ggplot2::labs(title = paste0("Density plot ", status, " data process") , subtitle = "")
    ## interactive
    for (i in seq(sample_list)) {
        if(i==1){
            density_sample <- abundance_long[which(
                abundance_long$sample==unique(abundance_long$sample)[i]),]
            density_data_x <- data.frame(density(density_sample$value)$x)
            density_data_y <- data.frame(density(density_sample$value)$y)
        }else{
            density_sample <- abundance_long[which(
                abundance_long$sample==unique(abundance_long$sample)[i]),]
            density_data_x <- data.frame(
                density_data_x,density(density_sample$value)$x)
            density_data_y <- data.frame(
                density_data_y,density(density_sample$value)$y)
        }
    }
    colnames(density_data_x) <- paste0("x", seq(sample_list))
    colnames(density_data_y) <- paste0("y", seq(sample_list))
    density_data <- cbind(density_data_x, density_data_y)
    for (i in seq(sample_list)) {
        if(i==1){
            densityPlot <- plotly::plot_ly(
                data = density_data, x=~x1, y=~y1, type = 'scatter',
                mode = 'lines', name = unique(abundance_long$sample)[1],
                hoverinfo="text", text=~paste(
                    "sample :", unique(abundance_long$sample)[1],
                    "\nabundance :", round(x1,3), "\ndensity :",
                    round(y1,3)))
        }else{
            text = paste0(
                "densityPlot %>% plotly::add_trace(x = ~x", i, ",", "y = ~y",i,
                ",name = unique(abundance_long$sample)[",i,
                "],text=~paste(\"sample :\",
                unique(abundance_long$sample)[",i,
                "],\"\nabundance :\",round(x",i,
                ",3),\"\ndensity :\",round(y",i,",3)))")
            densityPlot <- eval(parse(text=text))
        }
    }
    densityPlot <- densityPlot %>%
        plotly::layout(
            title = paste0("Density plot ", status, " data process"),
            legend=list(title=list(text="Sample")),
            xaxis=list(
                title = "Abundance",zeroline= FALSE,showgrid = FALSE),
            yaxis=list(title = "Density",showgrid = FALSE))
    return(list(densityPlot=densityPlot, density=density))
}

