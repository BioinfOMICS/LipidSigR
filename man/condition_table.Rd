\name{condition_table}
\alias{condition_table}
\docType{data}
\title{
condition_table
}
\description{
A data frame of sample names and clinical conditions for machine learning
analysis workflow. The data was modified from the paper "The landscape of cancer
cell line metabolism" and divided cancer cell lines into sensitive or resistant
to SCD gene knockout evenly based on gene dependency scores (CERES). 88 lipid
species were used to predict 228 cancer cell lines with a label of 0 (sensitive)
or 1 (resistant) by a binary classifier.
}
\usage{data("condition_table")}
\format{
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
data(ml_abundance)
}
\keyword{datasets}
