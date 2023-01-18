\name{DE_group_info}
\alias{DE_group_info}
\docType{data}
\title{
DE_group_info
}
\description{
A data frame of one clinical term (or one set of clinical terms) from
demographic data for differential expression analysis workflow. The grouping
information of samples. 10 healthy controls are considered as the control group
while 13 patients with systolic heart failure (HFrEF) are considered as
the experiment group.
}
\usage{data("DE_group_info")}
\format{
  A data frame with 23 observations on the following 4 variables.
  \describe{
    \item{\code{sample_name}}{a character vector, name of each sample.}
    \item{\code{label_name}}{a character vector, label name of each sample.}
    \item{\code{group}}{a character vector, group of each sample.}
    \item{\code{pair}}{a logical vector, paired sample or not.}
  }
}
\examples{
data(DE_group_info)
}
\keyword{datasets}
