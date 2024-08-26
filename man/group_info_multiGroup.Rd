\name{group_info_multiGroup}
\alias{group_info_multiGroup}
\docType{data}
\title{
group_info_multiGroup
}
\description{
A data frame of grouping information of samples
for multi-group differential expression analysis workflow.
This dataset includes three distinct
groups: a control group (N = 7) and two experimental groups: DHA-supplementation
(N = 6) and AA-supplementation (N = 4).
}
\usage{data("group_info_multiGroup")}
\format{
  A data frame with 17 observations on the following 3 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{label_name}}{a character vector}
    \item{\code{group}}{a character vector}
  }
}
\source{
Levental KR, Malmberg E, Symons JL, Fan YY, Chapkin RS, Ernst R, Levental I.
Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary
lipid perturbations to maintain cellular fitness. Nat Commun.
2020 Mar 12;11(1):1339.
}
\examples{
data(group_info_multiGroup)
}
\keyword{datasets}

