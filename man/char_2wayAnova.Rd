% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/char_2wayAnova.R
\name{char_2wayAnova}
\alias{char_2wayAnova}
\title{char_2wayAnova}
\usage{
char_2wayAnova(
  processed_se,
  ratio_transform = c("none", "log2"),
  char_transform = c("none", "log10", "square", "cube")
)
}
\arguments{
\item{processed_se}{A SummarizedExperiment object constructed by
\code{\link{as_summarized_experiment}} and processed by \code{\link{data_process}}.}

\item{ratio_transform}{Character. Method for transform the ratio-based abundance.
Allowed methods include "none" and "log2". Select 'none' to skip data transformation.
Default is \code{'log2'}.}

\item{char_transform}{Character. Method for transform the lipid
characteristics-based abundance. Allowed methods include "none", "log10",
"square", and "cube". Select 'none' to skip data transformation. Default is \code{'log10'}.}
}
\value{
Return 1 table with two-way ANOVA results.
}
\description{
This function calculates two-way ANOVA for all lipid characteristics.
}
\examples{
data("de_data_twoGroup")
processed_se_twoGroup <- data_process(
    se=de_data_twoGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
twoWayAnova_twoGroup <- char_2wayAnova(
    processed_se_twoGroup, ratio_transform='log2', char_transform='log10')

data("se_multiGroup")
processed_se_multiGroup <- data_process(
    se=se_multiGroup, exclude_missing=TRUE, exclude_missing_pct=70,
    replace_na_method='min', replace_na_method_ref=0.5,
    normalization='Percentage', transform='log10')
twoWayAnova_multiGroup <- char_2wayAnova(
    processed_se_multiGroup, ratio_transform='log2', char_transform='log10')
}
