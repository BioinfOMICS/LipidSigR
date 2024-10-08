% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-df.R
\name{as_summarized_experiment}
\alias{as_summarized_experiment}
\title{as_summarized_experiment}
\usage{
as_summarized_experiment(
  abundance,
  goslin_annotation,
  group_info = NULL,
  n_group = c("two", "multiple"),
  paired_sample = FALSE
)
}
\arguments{
\item{abundance}{A data frame of predictors includes features (such
as molecules or lipid classes) and their abundance in each sample. NAs are
allowed in all columns except the first. The first column must be named "feature"
and contain names of unique lipid species.}

\item{goslin_annotation}{A data frame of lipid characteristics from the
\code{\link[rgoslin]{parseLipidNames}} function of \pkg{\link{rgoslin}} package.
Please exclude lipids with the "Grammar" column marked as 'NOT_PARSEABLE'.}

\item{group_info}{A data frame containing group information. The required
columns are listed as below.
\enumerate{
\item sample_name: name of each sample. NAs are not allowed; the sample names
must be unique and match those in the abundance data.
\item label_name: label name of each sample. NAs are not allowed.
\item group: group name of each sample. NAs are not allowed.
\item pair: (only required for two-group data) If paired samples are contained,
they should be sequentially numbered from 1 to N without missing, blank, or
skipped numbers. If not, the values must all be marked as NA.
}}

\item{n_group}{Character. 'two' for inputting two-group data and
'multiple' for inputting multi-group data.}

\item{paired_sample}{Logical/NULL. For two-group data, enter TRUE/FALSE
to indicate if paired samples are included. For multi-group data, enter NULL.}
}
\value{
Return a SummarizedExperiment object.
}
\description{
This function construct a SummarizedExperiment object for analysis.
by input data frames.
}
\examples{
library(dplyr)
data("abundance_twoGroup")
data("group_info_twoGroup")
parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_twoGroup$feature)
recognized_lipid <- parse_lipid$Original.Name[which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
abundance <- abundance_twoGroup \%>\% dplyr::filter(feature \%in\% recognized_lipid)
goslin_annotation <- parse_lipid \%>\% dplyr::filter(Original.Name \%in\% recognized_lipid)
se <- as_summarized_experiment(abundance, goslin_annotation,
     group_info=group_info_twoGroup, n_group='two', paired_sample=FALSE)
}
