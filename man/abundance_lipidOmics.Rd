\name{abundance_lipidOmics}
\alias{abundance_lipidOmics}
\docType{data}
\title{
abundance_lipidOmics
}
\description{
A data frame of lipidomic abundance data for the multi-omics integration
workflow (\code{\link{as_MultiAssayExperiment}}. The lipid abundance data were
obtained from the lipidomic profiling of a hepatocellular carcinoma (HCC)
cohort of 98 patients recruited between 2011 and 2013. Samples are paired with
the matching transcriptomic profiling in \code{\link{abundance_rna}} and the
sample grouping in \code{\link{group_info_multiOmics}}.
}
\usage{data("abundance_lipidOmics")}
\format{
  A data frame with 333 observations (lipid features) on 60 variables.
  The first column, \code{feature}, is a character vector of lipid names.
  The remaining 59 columns are numeric vectors of lipid abundance values,
  one column per patient sample, named with de-identified sample IDs
  (e.g. \code{t_10126828}) matching \code{abundance_rna} and
  \code{group_info_multiOmics}.
}
\source{
Ma WL, et al. Ether-lipids accumulation promotes hepatocellular carcinoma
progression linked to PPAR\enc{α}{alpha} deficiency. J Biomed Sci.
2025;32(1):89. \doi{10.1186/s12929-025-01178-y}
}
\examples{
data(abundance_lipidOmics)
}
\keyword{datasets}
