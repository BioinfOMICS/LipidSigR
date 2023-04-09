\name{corr_adjusted_table}
\alias{corr_adjusted_table}
\docType{data}
\title{
corr_adjusted_table
}
\description{
A data frame with additional variables to adjust confounding effects for
correlation analysis workflow. Clinical terms for correcting in linear
regression or logistic regression analysis.
}
\usage{data("corr_adjusted_table")}
\value{
  A data frame with 129 observations on the following 6 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a numeric vector}
    \item{\code{Smoking}}{a numeric vector}
    \item{\code{BMI}}{a numeric vector}
    \item{\code{FEV1}}{a numeric vector}
  }
}
\references{
Bowler, Russell P., et al. "Plasma sphingolipids associated with chronic
obstructive pulmonary disease phenotypes." American journal of respiratory
and critical care medicine 191.3(2015):275-284.
}
\examples{
data(corr_adjusted_table)
}
\keyword{datasets}
