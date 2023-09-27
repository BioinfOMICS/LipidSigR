\name{DE_data}
\alias{DE_data}
\docType{data}
\title{
DE_data
}
\description{
The example data for differential expression analysis, which is the formal class
'SummarizedExperiment'. The data is obtained from the research of
adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left Ventricular failure. Human plasma lipidome from 10 healthy controls and 13
patients with systolic heart failure (HFrEF) were analyzed by MS-based shotgun
lipidomics.
}
\usage{data("DE_data")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, and rowData.
}
\value{
    \item{\code{assay}}{ A 208*23 matrix.}
    \item{\code{colData}}{ A data frame with 23 observations on the following 4
    variables.}
    \item{\code{rowData}}{ A data frame with 208 observations on the following
    10 variables.}
}
\source{
Salatzki J, Foryst-Ludwig A, Bentele K, Blumrich A, Smeir E, Ban Z, et al.
Adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left ventricular failure. PLoS genetics. 2018;14(1):e1007171.
}
\examples{
data(DE_data)
}
\keyword{datasets}
