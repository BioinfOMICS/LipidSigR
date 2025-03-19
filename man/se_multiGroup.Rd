\name{se_multiGroup}
\alias{se_multiGroup}
\docType{data}
\title{
se_multiGroup
}
\description{
The example data for differential expression analysis, which is the formal class
'SummarizedExperiment'. The lipid abundance data were obtained from the research of
Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary
lipid perturbations to maintain cellular fitness. This dataset includes three distinct
groups: a control group (N = 7) and two experimental groups: DHA-supplementation
(N = 6) and AA-supplementation (N = 4).
}
\usage{data("se_multiGroup")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, and rowData.
}
\value{
    \item{\code{assay}}{ A 1183*17 matrix.}
    \item{\code{colData}}{ A data frame with 17 observations on the following 3
    variables.}
    \item{\code{rowData}}{ A data frame with 1183 observations on the following
    72 variables.}
}
\source{
Levental KR, Malmberg E, Symons JL, Fan YY, Chapkin RS, Ernst R, Levental I.
Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary
lipid perturbations to maintain cellular fitness. Nat Commun.
2020 Mar 12;11(1):1339.
}
\examples{
data(se_multiGroup)
}
\keyword{datasets}
