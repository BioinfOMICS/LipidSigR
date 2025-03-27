\name{ml_se_sub}
\alias{ml_se_sub}
\docType{data}
\title{
ml_se_sub
}
\description{
The results of running the machine learning model function on subset example data,
which is the formal class 'SummarizedExperiment.'
}
\usage{data("ml_se_sub")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, rowData, and metadata.
}
\value{
    \item{\code{assay}}{ A 60*88 matrix.}
    \item{\code{colData}}{ A data frame with 60 observations on the following 2
    variables.}
    \item{\code{rowData}}{ A data frame with 88 observations on the following
    72 variables.}
    \item{\code{metadata}}{ A list with 15 objects, including machine learning
    parameters and results.}
}
\examples{
data(ml_se_sub)
}
\keyword{datasets}
