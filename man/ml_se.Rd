\name{ml_se}
\alias{ml_se}
\docType{data}
\title{
ml_se
}
\description{
The results of running the machine learning model function on example data,
which is the formal class 'SummarizedExperiment.'
}
\usage{data("ml_se")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, rowData, and metadata.
}
\value{
    \item{\code{assay}}{ A 228*88 matrix.}
    \item{\code{colData}}{ A data frame with 228 observations on the following 2
    variables.}
    \item{\code{rowData}}{ A data frame with 88 observations on the following
    72 variables.}
    \item{\code{metadata}}{ A list with 15 objects, including machine learning
    parameters and results.}
}
\examples{
data(ml_data)
}
\keyword{datasets}
