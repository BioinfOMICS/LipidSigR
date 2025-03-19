\name{abundance_multiGroup}
\alias{abundance_multiGroup}
\docType{data}
\title{
abundance_multiGroup
}
\description{
A data frame of lipid abundance data for multi-group differential abundance
analysis workflow. The lipid abundance data were obtained from the research of
Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary
lipid perturbations to maintain cellular fitness. This dataset includes three distinct
groups: a control group (N = 7) and two experimental groups: DHA-supplementation
(N = 6) and AA-supplementation (N = 4).
}
\usage{data("abundance_multiGroup")}
\format{
  A data frame with 2487 observations on the following 18 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{ctrl1}}{a numeric vector}
    \item{\code{ctrl2}}{a numeric vector}
    \item{\code{ctrl3}}{a numeric vector}
    \item{\code{ctrl4}}{a numeric vector}
    \item{\code{ctrl5}}{a numeric vector}
    \item{\code{ctrl6}}{a numeric vector}
    \item{\code{ctrl7}}{a numeric vector}
    \item{\code{DHA1}}{a numeric vector}
    \item{\code{DHA2}}{a numeric vector}
    \item{\code{DHA3}}{a numeric vector}
    \item{\code{DHA4}}{a numeric vector}
    \item{\code{DHA5}}{a numeric vector}
    \item{\code{DHA6}}{a numeric vector}
    \item{\code{AA1}}{a numeric vector}
    \item{\code{AA2}}{a numeric vector}
    \item{\code{AA3}}{a numeric vector}
    \item{\code{AA4}}{a numeric vector}
  }
}
\source{
Levental KR, Malmberg E, Symons JL, Fan YY, Chapkin RS, Ernst R, Levental I.
Lipidomic and biophysical homeostasis of mammalian membranes counteracts dietary
lipid perturbations to maintain cellular fitness. Nat Commun.
2020 Mar 12;11(1):1339.
}
\examples{
data(abundance_multiGroup)
}
\keyword{datasets}
