\name{corr_lipid_char_table}
\alias{corr_lipid_char_table}
\docType{data}
\title{
corr_lipid_char_table
}
\description{
A data frame of lipid characteristics for correlation analysis workflow. The
lipid characteristic data is obtained from the cohort interrogating the
associations of plasma sphingolipids with subphenotypes of COPD, including
airflow obstruction, emphysema, and frequent exacerbations.
}
\usage{data("corr_lipid_char_table")}
\value{
  A data frame with 69 observations on the following 2 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{class}}{a character vector}
  }
}
\source{
Bowler, Russell P., et al. "Plasma sphingolipids associated with chronic
obstructive pulmonary disease phenotypes." American journal of respiratory
and critical care medicine 191.3(2015):275-284.
}
\examples{
data(corr_lipid_char_table)
}
\keyword{datasets}
