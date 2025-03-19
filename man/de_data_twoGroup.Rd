\name{de_data_twoGroup}
\alias{de_data_twoGroup}
\docType{data}
\title{
de_data_twoGroup
}
\description{
The example data for differential expression analysis, which is the formal class
'SummarizedExperiment'. The data is obtained from the research of
adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left Ventricular failure. Human plasma lipidome from 10 healthy controls and 13
patients with systolic heart failure (HFrEF) were analyzed by MS-based shotgun
lipidomics.
}
\usage{data("de_data_twoGroup")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, and rowData.
}
\value{
    \item{\code{assay}}{ A 192*23 matrix.}
    \item{\code{colData}}{ A data frame with 23 observations on the following 4
    variables.}
    \item{\code{rowData}}{ A data frame with 192 observations on the following
    72 variables.}
}
\source{
Salatzki J, Foryst-Ludwig A, Bentele K, Blumrich A, Smeir E, Ban Z, et al.
Adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left ventricular failure. PLoS genetics. 2018;14(1):e1007171.
}
\examples{
data(de_data_twoGroup)
}
\keyword{datasets}
