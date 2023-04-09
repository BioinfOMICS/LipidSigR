\name{ML_condition_table}
\alias{ML_condition_table}
\docType{data}
\title{
ML_condition_table
}
\description{
A data frame of sample names and clinical conditions for machine learning
analysis workflow. A selected clinical term (or one set of clinical terms) from
demographic data with an assigned group number 0 or 1 to each sample.
}
\usage{data("ML_condition_table")}
\value{
  A data frame with 228 observations on the following 2 variables.
  \describe{
    \item{\code{sample_name}}{a character vector}
    \item{\code{group}}{a numeric vector}
  }
}
\source{
Li H, Ning S, Ghandi M, Kryukov GV, Gopal S, Deik A, et al. The landscape of
cancer cell line metabolism. Nature medicine. 2019;25(5):850.
}
\examples{
data(ML_condition_table)
}
\keyword{datasets}
