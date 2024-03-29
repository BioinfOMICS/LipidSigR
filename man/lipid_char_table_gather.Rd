% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipid_char_table_gather.R
\name{lipid_char_table_gather}
\alias{lipid_char_table_gather}
\title{lipid_char_table_gather}
\usage{
lipid_char_table_gather(exp_data_SE, char_var)
}
\arguments{
\item{exp_data_SE}{A SummarizedExperiment object contains information about 
various features, such as molecules, lipid class, etc., and their 
expression levels for each sample. The assay is a matrix representing the 
expression of lipid features across all samples. Missing values (NAs)  are 
allowed. The row data  corresponds to specific lipid features, such as class 
and total length. The first column's name must be "feature" (lipid species), 
and NAs are allowed for this data. The column data comprises sample names, 
sample labels, group names, and pair numbers that represent 'the pair' for 
conducting t-tests or Wilcoxon tests. NAs are allowed.}

\item{char_var}{A character string of the first lipid characteristic
selected by users from the column name of the rowData of \bold{exp_data_SE},
such as total length.}
}
\value{
Return a SummarizedExperiment object.
}
\description{
Select the user-selected characteristics, and make it become 
a long table.
}
\examples{
library(SummarizedExperiment)
data("profiling_data")
exp_data_SE <- profiling_data
lipid_char_table <- as.data.frame(SummarizedExperiment::rowData(exp_data_SE))
char_var <- colnames(lipid_char_table)[-1]
lipid_char_gather <- lipid_char_table_gather(exp_data_SE, 
    char_var = char_var[8])
}
