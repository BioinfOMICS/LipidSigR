\name{ML_lipid_char_table}
\alias{ML_lipid_char_table}
\docType{data}
\title{
ML_lipid_char_table
}
\description{
A data frame of lipid characteristics for machine learning analysis workflow.
The lipid characteristic data of 89 lipid species was modified from the paper
"The landscape of cancer cell line metabolism" and divided cancer cell lines
into sensitive or resistant to SCD gene knockout evenly based on gene dependency
scores (CERES).
}
\usage{data("ML_lipid_char_table")}
\value{
  A data frame with 89 observations on the following 6 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{class}}{a character vector}
    \item{\code{structural_category}}{a character vector}
    \item{\code{functional_category}}{a character vector}
    \item{\code{totallength}}{a numeric vector}
    \item{\code{totaldb}}{a numeric vector}
  }
}
\source{
Li H, Ning S, Ghandi M, Kryukov GV, Gopal S, Deik A, et al. The landscape of
cancer cell line metabolism. Nature medicine. 2019;25(5):850.
}
\examples{
data(ML_lipid_char_table)
}
\keyword{datasets}
