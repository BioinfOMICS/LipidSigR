\name{abundance_rna}
\alias{abundance_rna}
\docType{data}
\title{
abundance_rna
}
\description{
A data frame of transcriptomic (RNA) abundance data for the multi-omics
integration workflow (\code{\link{as_MultiAssayExperiment}}. The RNA
abundance data were obtained from the transcriptomic profiling of a
hepatocellular carcinoma (HCC) cohort of 98 patients recruited between 2011
and 2013. Samples are paired with the matching lipidomic profiling in
\code{\link{abundance_lipidOmics}} and the sample grouping in
\code{\link{group_info_multiOmics}}.
}
\usage{data("abundance_rna")}
\format{
  A data frame with 1118 observations (genes) on 60 variables.
  The first column, \code{gene_symbol}, is a character vector of gene
  symbols (renamed to \code{feature} internally by
  \code{\link{as_summarized_experiment}}). The remaining 59 columns are
  numeric vectors of RNA abundance values, one column per patient sample,
  named with de-identified sample IDs (e.g. \code{t_10126828}) matching
  \code{abundance_lipidOmics} and \code{group_info_multiOmics}.
}
\source{
Ma WL, et al. Ether-lipids accumulation promotes hepatocellular carcinoma
progression linked to PPAR\enc{α}{alpha} deficiency. J Biomed Sci.
2025;32(1):89. \doi{10.1186/s12929-025-01178-y}
}
\examples{
data(abundance_rna)
}
\keyword{datasets}
