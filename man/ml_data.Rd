\name{ml_data}
\alias{ml_data}
\docType{data}
\title{
ml_data
}
\description{
The example data for machine learning analysis, which is the formal class
'SummarizedExperiment'. The data is obtained from the paper "The landscape of cancer
cell line metabolism" and divided cancer cell lines into sensitive or resistant
to SCD gene knockout evenly based on gene dependency scores (CERES). 88 lipid
species were used to predict 228 cancer cell lines with a label of 0 (sensitive)
or 1 (resistant) by a binary classifier.
}
\usage{data("ml_data")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, and rowData.
}
\value{
    \item{\code{assay}}{ A 228*88 matrix.}
    \item{\code{colData}}{ A data frame with 228 observations on the following 2
    variables.}
    \item{\code{rowData}}{ A data frame with 88 observations on the following
    72 variables.}
}
\source{
Li H, Ning S, Ghandi M, Kryukov GV, Gopal S, Deik A, et al. The landscape of
cancer cell line metabolism. Nature medicine. 2019;25(5):850.
}
\examples{
data(ml_data)
}
\keyword{datasets}
