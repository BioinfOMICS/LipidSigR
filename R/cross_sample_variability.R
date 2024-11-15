#' @title cross_sample_variability
#' @description This function provides three distribution plots for users to view sample
#' variability and compare lipids' amount/abundance difference between samples (i.e., patients vs. control).
#' @param se A SummarizedExperiment object construct by \code{\link{as_summarized_experiment}}.
#' @return Return a list of 3 interactive plots, 3 static plots, and 2 tables.
#' \enumerate{
#' \item interactive_lipid_number_barPlot & static_lipid_number_barPlot: histogram of lipid numbers.
#' \item interactive_lipid_amount_barPlot & static_lipid_amount_barPlot: histogram of the total amount of lipid in each sample.
#' \item interactive_lipid_distribution & static_lipid_distribution: density plot of the underlying probability distribution of the lipid abundance in each sample.
#' \item table_total_lipid: table for plotting lipid_number_barPlot & lipid_amount_barPlot.
#' \item table_lipid_distribution: table for plotting lipid_distribution.
#' }
#' @export
#' @examples
#' data("profiling_data")
#' result <- cross_sample_variability(profiling_data)

cross_sample_variability <- function(se){
    # check input SE
    .check_inputSE(se, metadata_list=NULL)

    abundance <- .extract_df(se, type = "abundance")
    abundance_trans_data <- abundance %>%
        tidyr::gather(sample_name, value, -feature)
    num.lip <- abundance_trans_data %>%
        dplyr::mutate(is.abundance=!is.na(value)) %>%
        dplyr::group_by(sample_name) %>%
        dplyr::summarise(abundance.count=sum(is.abundance)) %>%
        dplyr::arrange(sample_name)
    num.lip$sample_name <- factor(
        num.lip$sample_name, levels=unique(num.lip$sample_name)[order(
            num.lip$abundance.count, decreasing=TRUE)])
    tot.lip <- abundance_trans_data %>%
        dplyr::group_by(sample_name) %>%
        dplyr::summarise(lipid_amount=sum(value, na.rm=TRUE))
    tot.lip$sample_name <- factor(
        tot.lip$sample_name, levels=unique(tot.lip$sample_name)[order(
            tot.lip$lipid_amount, decreasing=TRUE)])
    tot.num.lip <- merge(x=num.lip, y=tot.lip, by="sample_name", all=TRUE)
    dens.lip <- stats::na.omit(abundance_trans_data)

    ## plotting
    dynamic_plot <- .dynamic_sample_variability(tot.num.lip, dens.lip)
    static_plot <- .static_sample_variability(tot.num.lip, dens.lip)

    return(list(
        interactive_lipid_number_barPlot = dynamic_plot$abundance_lip,
        interactive_lipid_amount_barPlot = dynamic_plot$lipid_amount,
        interactive_lipid_distribution = dynamic_plot$abundance_density,
        static_lipid_number_barPlot = static_plot$abundance_lip,
        static_lipid_amount_barPlot = static_plot$lipid_amount,
        static_lipid_distribution = static_plot$abundance_density,
        table_total_lipid=tot.num.lip, table_lipid_distribution=dens.lip))
}

## plotly total lipid amount of samples & hist of quantification per sample
.dynamic_sample_variability <- function(tot.num.lip, dens.lip){
    abundance_lip <- plotly::plot_ly(
        tot.num.lip, x=~sample_name,y=~abundance.count, type="bar",
        color = ~sample_name,hoverinfo="text", text=~paste(
            "Sample name :",sample_name,"\nCount :",abundance.count)) %>%
        plotly::layout(legend = list(
            title=list(text="Sample"), font=list(family='arial')),
            xaxis = list(
                titlefont=list(family='arial'), tickfont=list(family='arial'),
                automargin = TRUE, showgrid=FALSE),
            yaxis = list(
                title="Number of Expressed Lipids",
                titlefont=list(family='arial'), tickfont=list(family='arial'),
                automargin=TRUE, showgrid=FALSE))
    lipid_amount <- plotly::plot_ly(
        tot.num.lip, x=~sample_name,y=~lipid_amount, type="bar",
        color=~sample_name, hoverinfo="text",
        text=~paste(
            "Sample name :", sample_name, "\nAmount :",round(lipid_amount,3))) %>%
        plotly::layout(legend = list(
            title=list(text="Sample"), font = list(family = 'arial')),
            xaxis = list(
                titlefont=list(family='arial'), tickfont=list(family='arial'),
                automargin=TRUE, showgrid=FALSE),
            yaxis = list(
                title="Lipid Amount", titlefont=list(family='arial'),
                tickfont=list(family='arial'), automargin=TRUE, showgrid=FALSE))
    # density
    length_sample <- length(unique(dens.lip$sample_name))
    for (i in 1:length_sample) {
        if(i==1){
            density_sample <- dens.lip[which(
                dens.lip$sample_name==unique(dens.lip$sample_name)[i]),]
            density_data_x <- data.frame(density(log10(density_sample$value))$x)
            density_data_y <- data.frame(density(log10(density_sample$value))$y)
        }else{
            density_sample <- dens.lip[which(
                dens.lip$sample_name==unique(dens.lip$sample_name)[i]),]
            density_data_x <- data.frame(
                density_data_x,density(log10(density_sample$value))$x)
            density_data_y <- data.frame(
                density_data_y,density(log10(density_sample$value))$y)
        }
    }
    colnames(density_data_x) <- paste0("x",1:length_sample)
    colnames(density_data_y) <- paste0("y",1:length_sample)
    density_data <- cbind(density_data_x,density_data_y)
    for (i in 1:length_sample) {
        if(i==1){
            abundance_density <- plotly::plot_ly(
                data=density_data, x=~x1, y=~y1, type='scatter',
                mode='lines', name=unique(dens.lip$sample_name)[1], hoverinfo="text",
                text=~paste(
                    "Sample name :", unique(dens.lip$sample_name)[1],
                    "\nlog10(abundance) :", round(x1,3), "\nDensity :", round(y1,3)))
        }else{
            text = paste0(
                "abundance_density %>% plotly::add_trace(x = ~x", i, ",", "y = ~y",i,
                ",name = unique(dens.lip$sample_name)[",i, "],text=~paste(\"Sample name :\",
                unique(dens.lip$sample_name)[",i, "],\"\nlog10(abundance) :\",round(x",i,
                ",3),\"\nDensity :\",round(y",i,",3)))")
            abundance_density <- eval(parse(text=text))
        }
    }
    abundance_density <- abundance_density %>%
        plotly::layout(
            legend=list(title=list(text="Sample")),
            xaxis=list(
                title="log10(abundance)", zeroline= FALSE, showgrid=FALSE),
            yaxis=list(title="Density", showgrid=FALSE))
    return(list(abundance_lip=abundance_lip, lipid_amount=lipid_amount,
                abundance_density=abundance_density))
}

.static_sample_variability <- function(tot.num.lip, dens.lip){
    abundance_lip <- ggplot2::ggplot(
        tot.num.lip, ggplot2::aes(
            x=sample_name, y=abundance.count, fill=sample_name)) +
        ggplot2::geom_col() +
        ggplot2::labs(fill="Sample", y="Number of Expressed Lipids", x="") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle=90, vjust=0.5,hjust=1)) +
        ggthemes::theme_hc()
    lipid_amount <- ggplot2::ggplot(
        tot.num.lip, ggplot2::aes(
            x=sample_name, y=lipid_amount, fill=sample_name)) +
        ggplot2::geom_col() +
        ggplot2::labs(fill='Sample', y='Lipid Amount', x="") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1)) +
        ggthemes::theme_hc()
    abundance_density <- ggplot2::ggplot(
        dens.lip, ggplot2::aes(x=log10(value), color=sample_name)) +
        ggplot2::geom_density() +
        ggplot2::labs(color='Sample', y='Density', x="log10(abundance)") +
        ggthemes::theme_hc()
    return(list(abundance_lip=abundance_lip, lipid_amount=lipid_amount,
                abundance_density=abundance_density))
}
