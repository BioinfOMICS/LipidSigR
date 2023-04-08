\name{DE_lipid_char_table}
\alias{DE_lipid_char_table}
\docType{data}
\title{
DE_lipid_char_table
}
\description{
A data frame of lipid characteristics for differential expression analysis
workflow. The lipid characteristic data were obtained from the research of
adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left Ventricular failure.
}
\usage{data("DE_lipid_char_table")}
\value{
  A data frame with 208 observations on the following 10 variables.
  \describe{
    \item{\code{feature}}{a character vector}
    \item{\code{class}}{a character vector}
    \item{\code{structural_category}}{a character vector}
    \item{\code{functional_category}}{a character vector}
    \item{\code{totallength}}{a numeric vector}
    \item{\code{totaldb}}{a numeric vector}
    \item{\code{totaloh}}{a numeric vector}
    \item{\code{FA_length}}{a numeric vector}
    \item{\code{FA_db}}{a character vector}
    \item{\code{FA_oh}}{a character vector}
  }
}
\source{
Salatzki J, Foryst-Ludwig A, Bentele K, Blumrich A, Smeir E, Ban Z, et al.
Adipose tissue ATGL modifies the cardiac lipidome in pressure-overload-induced
left ventricular failure. PLoS genetics. 2018;14(1):e1007171.
}
\examples{
data(DE_lipid_char_table)
}
\keyword{datasets}
