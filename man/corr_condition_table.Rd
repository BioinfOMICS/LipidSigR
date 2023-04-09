\name{corr_condition_table}
\alias{corr_condition_table}
\docType{data}
\title{
corr_condition_table
}
\description{
A data frame of sample names and clinical conditions for correlation analysis
workflow. A table of desired clinical terms (i.e., clinical subphenotypes of
chronic obstructive pulmonary disease).
}
\usage{data("corr_condition_table")}
\value{
  A data frame with 129 observations on the following 4 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{FEV1_FVC}}{a numeric vector}
    \item{\code{Emphysema}}{a numeric vector}
    \item{\code{Exacerbations}}{a numeric vector}
  }
}
\source{
Bowler, Russell P., et al. "Plasma sphingolipids associated with chronic
obstructive pulmonary disease phenotypes." American journal of respiratory and
critical care medicine 191.3(2015):275-284.
}
\examples{
data(corr_condition_table)
}
\keyword{datasets}
