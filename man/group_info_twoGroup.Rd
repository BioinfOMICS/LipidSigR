\name{group_info_twoGroup}
\alias{group_info_twoGroup}
\docType{data}
\title{
group_info_twoGroup
}
\description{
A data frame of grouping information of samples
for two-group differential expression analysis workflow. 10 healthy controls
are considered as the control group while 13 patients with systolic
heart failure (HFrEF) are considered as the experiment group.
}
\usage{data("group_info_twoGroup")}
\format{
  A data frame with 23 observations on the following 4 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{label_name}}{a character vector}
    \item{\code{group}}{a character vector}
    \item{\code{pair}}{a logical vector}
  }
}
\source{
Salatzki J, Foryst-Ludwig A, Bentele K, Blumrich A, Smeir E, Ban Z, et al.
Adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left ventricular failure. PLoS genetics. 2018;14(1):e1007171.
}
\examples{
data(group_info_twoGroup)
}
\keyword{datasets}
