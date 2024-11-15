#' @title lipid_profiling
#' @description This function provides two plots for users to explore lipid
#' abundance over a specific lipid characteristic.
#' @param processed_se A SummarizedExperiment object constructed by
#' \code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.
#' @param char Character. A lipid characteristic selected from the common list r
#' eturned by \code{\link{list_lipid_char}}.
#' @return Return a list of 2 interactive plots, 2 static plots, and 2 tables.
#' \enumerate{
#' \item interactive_char_barPlot & static_char_barPlot: bar plot classified by
#' selected lipid characteristic.
#' within each group (e.g., PE, PC) of selected characteristics (e.g., class).
#' \item interactive_lipid_composition & static_lipid_composition: stacked horizontal
#' bar chart of lipid class composition.
#' \item table_char_barPlot: table for plotting char_barPlot.
#' \item table_lipid_composition: table for plotting lipid_composition.
#' }
#' @export
#' @examples
#' data("profiling_data")
#' processed_se <- data_process(
#'     profiling_data, exclude_missing=TRUE, exclude_missing_pct=70,
#'     replace_na_method='min', replace_na_method_ref=0.5,
#'     normalization='Percentage')
#' list_lipid_char(processed_se)$common_list
#' result <- lipid_profiling(processed_se, char="class")

lipid_profiling <- function(processed_se, char){
    # check input SE
    .check_inputSE(processed_se, metadata_list=NULL)
    if (is.null(char) | isFALSE(.check_char(processed_se, char, type='common'))) {
        stop("Wrong char input, you can view the available char list by list_lipid_char function.")
    }
    abundance <- .extract_df(processed_se, type = "abundance")
    .check_imputation(abundance)

    ## lipid composition
    lipid_char_table_raw <- .extract_df(processed_se, type = "lipid")
    char_list <- list_lipid_char(processed_se)$common_list
    names(char_list) <- NULL
    char_col <- c(TRUE,
                  colnames(lipid_char_table_raw)[-1] %in% char_list)
    lipid_char_table <- lipid_char_table_raw[, char_col]
    lipid_char_gather <- lipid_char_table %>%
        tidyr::gather(lipid.category, lipid.category.value, -feature)
    ## process "|" characteristic
    if(any(grepl("\\|", lipid_char_gather)  )){
        lipid_char_table_split <- tidyr::separate_rows(
            lipid_char_gather, "lipid.category.value", sep = "\\|")
    } else {
        lipid_char_table_split <- lipid_char_gather
    }

    lipid_composion <- merge(
        abundance %>% tidyr::gather(sample_name, value, -feature),
        lipid_char_table_split, by="feature") %>%
        dplyr::filter(
            lipid.category==char & !is.na(lipid.category.value)) %>%
        dplyr::group_by(
            sample_name, lipid.category, lipid.category.value) %>%
        dplyr::summarise(value=sum(value, na.rm=TRUE)) %>%
        plotly::ungroup() %>% dplyr::group_by(sample_name) %>%
        dplyr::mutate(weight=100/sum(value)) %>%
        dplyr::mutate(value=value*weight) %>% plotly::ungroup()

    dynamic_plot <- .dynamic_profiling(lipid_composion, char)
    static_plot <- .static_profiling(lipid_composion, char)

    return(list(
        interactive_char_barPlot=dynamic_plot$barPlot,
        interactive_lipid_composition=dynamic_plot$composion_plot,
        static_char_barPlot=static_plot$barPlot,
        static_lipid_composition=static_plot$composion_plot,
        table_char_barPlot=lipid_composion,
        table_lipid_composition=lipid_composion))
}

.dynamic_profiling <- function(lipid_composion, char){
    barPlot <- lipid_composion %>%
        plotly::plot_ly(
            x=~lipid.category.value, y=~value, color = ~sample_name, type="bar",
            hoverinfo="text", text=~paste(
                Hmisc::capitalize(char)," :",lipid.category.value,
                "\nSample name :",sample_name, "\nLipid Abundance :",
                round(value,3))) %>%
        plotly::layout(legend = list(title=list(text="Sample name")),
               xaxis = list(title=Hmisc::capitalize(char)),
               yaxis = list(title='Lipid Abundance',nticks=15),
               title = paste0('Lipid ', Hmisc::capitalize(char)))
    composion_plot <- lipid_composion %>%
        plotly::plot_ly(
            x=~value, y=~sample_name, color = ~lipid.category.value, type="bar",
            hoverinfo="text", orientation = "h", text=~paste(
                Hmisc::capitalize(char),":",lipid.category.value,
                "\nSample name :",sample_name, "\nLipid Abundance :",
                round(value,3))) %>%
        plotly::layout(
            barmode = 'stack', legend = list(
                title=list(text=Hmisc::capitalize(char))),
            xaxis = list(title="%",nticks=15),
            yaxis = list(title=''),
            title = paste0(
                'Lipid ', Hmisc::capitalize(char), ' Composition'))
    return(list(barPlot=barPlot, composion_plot=composion_plot))
}

.static_profiling <- function(lipid_composion, char){
    barPlot <- ggplot2::ggplot(
        lipid_composion, ggplot2::aes(
            x=lipid.category.value, y=value, fill=sample_name)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::labs(
            fill="Sample name", y='Lipid Abundance',
            x=Hmisc::capitalize(char),
            title=paste0('Lipid ', Hmisc::capitalize(char))) +
        ggthemes::theme_hc()
    composion_plot <- ggplot2::ggplot(
        lipid_composion, ggplot2::aes(
            x=value, y=sample_name, fill=lipid.category.value)) +
        ggplot2::geom_col() +
        ggplot2::labs(
            fill=Hmisc::capitalize(char), y='', x="%",
            title=paste0('Lipid ',
                         Hmisc::capitalize(char), 'Composition')) +
        ggthemes::theme_hc()
    return(list(barPlot=barPlot, composion_plot=composion_plot))
}
