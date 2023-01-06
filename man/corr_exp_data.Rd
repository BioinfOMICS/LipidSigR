\name{corr_exp_data}
\alias{corr_exp_data}
\docType{data}
\title{
corr_exp_data
}
\description{
A data frame of lipid expression data for correlation analysis workflow.
}
\usage{data("corr_exp_data")}
\format{
  A data frame with 69 observations on the following 130 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{sample1-129}}{a numeric vector, expression value}
  }
}
\details{
The lipid expression data is obtained from the cohort interrogating the
associations of plasma sphingolipids with subphenotypes of COPD,
including airflow obstruction, emphysema, and frequent exacerbations.
The dataset contained 129 current and former smokers.
}
\source{
Bowler, Russell P., et al. "Plasma sphingolipids associated with chronic
obstructive pulmonary disease phenotypes." American journal of respiratory
and critical care medicine 191.3 (2015): 275-284.
}
\examples{
data(corr_exp_data)
}
\keyword{datasets}
