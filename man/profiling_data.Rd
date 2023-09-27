\name{profiling_data}
\alias{profiling_data}
\docType{data}
\title{
profiling_data
}
\description{
The example data for profiling analysis, which is the formal class
'SummarizedExperiment'. The data is obtained from the research of adipose
tissue ATGL modifies the cardiac lipidome in pressure-overload-induced left
Ventricular failure. Human plasma lipidome from 10 healthy controls and 13
patients with systolic heart failure (HFrEF) were analyzed by MS-based shotgun
lipidomics. The lipid characteristic data is obtained from the cohort
interrogating the associations of plasma sphingolipids with subphenotypes of
COPD, including airflow obstruction, emphysema, and frequent exacerbations.
}
\usage{data("profiling_data")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, rowData, and colData.
}
\value{
    \item{\code{assay}}{ A 208*23 matrix.}
    \item{\code{rowData}}{ A data frame with 208 observations on the following
    10 variables.}
    \item{\code{colData}}{ A data frame with 23 observations on the following 1
    variables.}
}
\source{
Salatzki J, Foryst-Ludwig A, Bentele K, Blumrich A, Smeir E, Ban Z, et al.
Adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left ventricular failure. PLoS genetics. 2018;14(1):e1007171.
}
\examples{
data(profiling_data)
}
\keyword{datasets}
