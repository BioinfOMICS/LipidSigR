\name{corr_data}
\alias{corr_data}
\docType{data}
\title{
corr_data
}
\description{
The example data for correlation analysis, which is the formal class
'SummarizedExperiment'.
The data was obtained from Bowler RP et al.'s study, "Plasma Sphingolipids
Associated with Chronic Obstructive Pulmonary Disease Phenotypes." Targeted mass
spectrometry identified 69 distinct plasma sphingolipid species in 129 current
and former smokers. This cohort was analyzed to explore the associations between
plasma sphingolipids and COPD subphenotypes, including airflow obstruction,
emphysema, and frequent exacerbations.
}
\usage{data("corr_data")}
\format{
The format is: Formal class 'SummarizedExperiment'
[package "SummarizedExperiment"] with assay, colData, and rowData.
}
\value{
    \item{\code{assay}}{ A 50*129 matrix.}
    \item{\code{colData}}{ A data frame with 129 observations on the following 9
    variables.}
    \item{\code{rowData}}{ A data frame with 50 observations on the following
    72 variables.}
}
\source{
Bowler RP, Jacobson S, Cruickshank C, Hughes GJ, Siska C, Ory DS, Petrache I, Schaffer JE, Reisdorph N, Kechris K. Plasma sphingolipids associated with chronic obstructive pulmonary disease phenotypes. Am J Respir Crit Care Med. 2015;191(3):275-84.
}
\examples{
data(corr_data)
}
\keyword{datasets}
